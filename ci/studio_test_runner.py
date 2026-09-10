"""Report desktop failures immediately and bound each test, including Qt teardown."""

import faulthandler
import traceback
import unittest


class DesktopResult(unittest.TextTestResult):
    def startTest(self, test):
        # The suite grows as workspaces are added; a progressing suite must not
        # consume a single two-minute allowance shared by every test.
        faulthandler.dump_traceback_later(120, exit=True)
        super().startTest(test)

    def addError(self, test, error):
        super().addError(test, error)
        self._report_now(test, error)

    def addFailure(self, test, error):
        super().addFailure(test, error)
        self._report_now(test, error)

    def _report_now(self, test, error):
        # A later native crash must not hide failures already detected.
        self.stream.writeln(f"\nImmediate failure details: {test.id()}")
        self.stream.write("".join(traceback.format_exception(*error)))
        self.stream.flush()


class DesktopRunner(unittest.TextTestRunner):
    resultclass = DesktopResult


if __name__ == "__main__":
    faulthandler.enable()
    faulthandler.dump_traceback_later(120, exit=True)  # Includes test discovery.
    try:
        unittest.main(module=None, testRunner=DesktopRunner)
    finally:
        faulthandler.cancel_dump_traceback_later()
