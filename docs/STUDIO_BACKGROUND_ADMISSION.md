# Background CAD and archive admission

Studio `0.6.0a2.dev5` separates response reading and numerical validation from
the GUI event loop. This extends the [shared CPU admission](ENGINE_CPU_ADMISSION.md)
introduced in dev4; the native kernel remains version 0.20.0.

| Operation | Background work | Publication on the GUI thread |
|---|---|---|
| CAD creation or feature preview | Bounded JSON read, body decoding, geometry/mass/frame/recipe and CPU checks, immutable document preparation with updated attachments | Adopt the prepared document, commit undo history, update the scene |
| Open CalculiX study or result | Existing study/result loader and its numerical consistency checks | Show the accepted study or archived result |
| Open articulated input or result | Existing input/operator contract, units, state and numerical consistency checks | Adopt an accepted input after the discard prompt, or show the archived operators |

CAD admission retains the operation's CPU lease until its reader exits. Stored
archive reads request one slot from the same FIFO scheduler; a queued read does
not touch its files. The reader receives plain data and an immutable project
snapshot. Qt widgets and VTK scene mutation stay on the GUI thread, with bound
receiver slots handling completion.

Cancel requests interruption and suppresses publication. A reader already inside
a parsing or numerical routine may finish that routine before stopping. Closing
a window or CAD dialog defers its destruction until then; it does not forcibly
terminate a thread. Temporary CAD files and the CPU allocation remain available
to the reader. Previous accepted results and pending analysis inputs survive
failure or cancellation. CAD publication also checks that the captured document
is still the current document. Undo and redo preserve their fresh revisions.

The [real Qt tests](../apps/studio/tests/test_background_admission.py) cover a
queued read without file access, cancellation at launch, an event-loop heartbeat
while a real result loader is held at a controlled barrier, preserved edited
inputs and displayed results, deferred destruction in both analysis windows and
both CAD dialogs, and a real CAD command that rejects body decoding on the GUI
thread. Existing CAD rejection, frame/attachment, archive and workspace tests
remain applicable.

The [retained installed-wheel qualification](bancs/background-admission-dev5/README.md)
includes the complete Studio suite, separate Pinocchio tests, keyboard recipes,
normal startup and ten repetitions of the new lifecycle tests in one Qt process.

These are lifecycle and thread-ownership checks, not latency or throughput
guarantees for arbitrary files. Python validation still shares the interpreter
GIL. CAD request serialization, ordinary project file I/O, STEP export writing,
history revision validation and VTK scene updates retain GUI work. This change
does not claim that every expensive operation has moved off the GUI thread or
that the scheduler limits every numerical-library or operating-system thread.
