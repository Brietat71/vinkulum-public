"""Cancellable archive reads with CPU admission and receiver-bound Qt lifecycle."""

from PySide6.QtCore import QObject, QThread, Signal, Slot

from .cpu_scheduler import application_scheduler


class ArchiveRead(QThread):
    def __init__(self, loader, arguments, context, parent):
        super().__init__(parent)
        self.loader, self.arguments, self.context = loader, arguments, context
        self.value, self.error = None, None

    def run(self):
        try:
            value = self.loader(*self.arguments)
            if not self.isInterruptionRequested():
                self.value = value
        except (OSError, ValueError, TypeError, KeyError, RuntimeError) as error:
            self.error = str(error)


class ArchiveController(QObject):
    busy_changed = Signal(bool)
    stage_changed = Signal(str)
    completed = Signal(object, object)
    problem = Signal(str)
    cancelled = Signal(str)

    def __init__(self, parent=None, *, scheduler=None):
        super().__init__(parent)
        self.scheduler = scheduler if scheduler is not None else application_scheduler()
        self.busy, self.reader, self._ticket, self._reason = False, None, None, None

    def start(self, loader, *arguments, context=None):
        if self.busy:
            raise RuntimeError("An archive is already being checked.")
        if not callable(loader):
            raise TypeError("An archive reader is required.")
        self._job = (loader, arguments, context)
        self._reason = None
        self.busy = True
        ticket = self.scheduler.request(1)
        self._ticket = ticket
        ticket.ready.connect(self._allocated)
        self.busy_changed.emit(True)
        if ticket is self._ticket:
            self.stage_changed.emit("Waiting for archive validation CPU resources…")

    @Slot(object)
    def _allocated(self, plan):
        if not self.busy or self._reason:
            return
        reader = ArchiveRead(*self._job, self)
        self.reader = reader
        reader.finished.connect(self._reader_finished)
        self.stage_changed.emit(
            "Reading archive and checking its numerical consistency…"
        )
        if self.reader is not reader:
            return
        if self._reason:
            self._finished(reader)
        else:
            reader.start()

    def cancel(self, reason="Archive loading cancelled. Previous data preserved."):
        if not self.busy:
            return
        self._reason = reason
        if self.reader is not None:
            self.reader.requestInterruption()
        else:
            self._settle(reason, cancelled=True)

    @Slot()
    def _reader_finished(self):
        self._finished(self.sender())

    def _finished(self, reader):
        if reader is None or reader is not self.reader:
            return
        self.reader = None
        reader.deleteLater()
        reason = self._reason or reader.error
        if reason or reader.value is None:
            self._settle(
                reason or "Archive returned no validated data.",
                cancelled=self._reason is not None,
            )
            return
        value, context = reader.value, reader.context
        self._settle()
        self.completed.emit(value, context)

    def _settle(self, reason=None, *, cancelled=False):
        if not self.busy:
            return
        self.busy = False
        self._job = None
        if self._ticket is not None:
            self._ticket.release()
            self._ticket = None
        self.busy_changed.emit(False)
        if reason:
            (self.cancelled if cancelled else self.problem).emit(reason)

    def shutdown(self):
        self.cancel("Window closed. Archive loading stopped.")
        reader = self.reader
        if reader is not None:
            if not reader.wait(2000):
                raise RuntimeError(
                    "Archive checks are still stopping; keep their owner alive."
                )
            self._finished(reader)
