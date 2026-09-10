"""Qt-thread CPU admission shared by native and external Studio engines.

Tokens bound admitted compute allocations, not every operating-system thread.
Dispatch is deferred so signal handlers can cancel a job before its process starts.
"""

from PySide6.QtCore import QCoreApplication, QObject, QTimer, Signal

from .execution import ExecutionPlan, default_budget


class CpuTicket(QObject):
    ready = Signal(object)

    def __init__(self, scheduler, threads):
        super().__init__(scheduler)
        self.scheduler = scheduler
        self.plan = ExecutionPlan(threads, scheduler.capacity)
        self.state = "queued"

    def release(self):
        self.scheduler.release(self)


class CpuScheduler(QObject):
    def __init__(self, capacity, parent=None):
        super().__init__(parent if parent is not None else QCoreApplication.instance())
        ExecutionPlan(capacity, capacity)
        self.capacity = capacity
        self._queue = []
        self._active = set()

    @property
    def allocated(self):
        return sum(ticket.plan.threads for ticket in self._active)

    def request(self, threads=None):
        ticket = CpuTicket(self, self.capacity if threads is None else threads)
        self._queue.append(ticket)
        QTimer.singleShot(0, self._dispatch)
        return ticket

    def release(self, ticket):
        if ticket.state == "released":
            return
        if ticket in self._queue:
            self._queue.remove(ticket)
        self._active.discard(ticket)
        ticket.state = "released"
        ticket.deleteLater()
        QTimer.singleShot(0, self._dispatch)

    def _dispatch(self):
        while self._queue:
            ticket = self._queue[0]
            if self.allocated + ticket.plan.threads > self.capacity:
                return
            self._queue.pop(0)
            self._active.add(ticket)
            ticket.state = "running"
            ticket.ready.emit(ticket.plan)


def application_scheduler():
    app = QCoreApplication.instance()
    if app is None:
        raise RuntimeError("Studio CPU scheduling requires a Qt application.")
    if not hasattr(app, "_vinkulum_cpu_scheduler"):
        app._vinkulum_cpu_scheduler = CpuScheduler(default_budget(), app)
    return app._vinkulum_cpu_scheduler
