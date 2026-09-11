"""Build the qualified FreeCAD 1.1.3 host in a new Linux x86-64 directory.

Requires Git and Pixi 0.80.0. No system installation is performed.
SPDX-License-Identifier: Apache-2.0
"""

import argparse
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
from contextlib import contextmanager
from pathlib import Path

COMMIT = '145529fe741292ff0b3977a01195bf0247425794'
LOCK = 'f9a9c1f90cd7aa37655c593793c67dacec3e7b532b2403a29ff2575f013947f9'
PATCHES = ('freecad-1.1.3-brep-floatfield.patch', 'freecad-1.1.3-placement-roundtrip.patch')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@contextmanager
def destination_lock(root):
    import fcntl

    descriptor = os.open(root / '.build.lock', os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW, 0o600)
    try:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ValueError('Another installer owns this destination; wait for it to finish.') from error
        yield
    finally:
        os.close(descriptor)


def require_idle(report):
    for step in report.get('steps', []):
        if step.get('status') == 'passed':
            continue
        group = step.get('process_session')
        if group is None:
            raise ValueError('An interrupted launch has no process record; use a new destination.')
        # Ninja may create separate process groups for commands in the same session.
        for proc in Path('/proc').iterdir():
            if not proc.name.isdigit():
                continue
            try:
                fields = (proc / 'stat').read_text().rsplit(')', 1)[1].split()
            except FileNotFoundError:
                continue
            if int(fields[3]) == group and fields[0] != 'Z':
                raise ValueError('A previous build session is still active; wait for it to finish.')


def run_recorded(command, root, environment, log, step, save):
    step['status'] = 'starting'
    save()
    process = subprocess.Popen(command, cwd=root, env=environment, stdout=log,
                               stderr=subprocess.STDOUT, start_new_session=True)
    step.update(status='running', process_session=process.pid)
    save()
    code = process.wait()
    require_idle({'steps': [step]})
    if code:
        raise subprocess.CalledProcessError(code, command)
    step['status'] = 'passed'
    save()


def build(args):
    if platform.system() != 'Linux' or platform.machine() != 'x86_64':
        raise ValueError('This host recipe is qualified only for Linux x86-64.')
    if not 1 <= args.jobs <= 32:
        raise ValueError('Choose 1–32 build jobs.')
    root = args.directory.expanduser().absolute()
    if root.is_symlink():
        raise ValueError('The destination must not be a symbolic link.')
    pixi = shutil.which(args.pixi)
    if pixi is None or shutil.which('git') is None:
        raise ValueError('Git and Pixi 0.80.0 are required.')
    if subprocess.check_output([pixi, '--version'], text=True).strip() != 'pixi 0.80.0':
        raise ValueError('Use the qualified Pixi 0.80.0 version.')
    patches = Path(__file__).resolve().parent / 'patches'
    expected = {name: digest(patches / name) for name in PATCHES}
    environment = os.environ.copy()
    for name in ('PYTHONHOME', 'PYTHONPATH', 'VIRTUAL_ENV', 'LD_LIBRARY_PATH',
                 'LD_PRELOAD', 'QT_PLUGIN_PATH', 'QT_QPA_PLATFORM_PLUGIN_PATH'):
        environment.pop(name, None)
    environment.update(PIXI_HOME=str(root / 'pixi-home'), PIXI_CACHE_DIR=str(root / 'cache'),
                       CCACHE_DIR=str(root / 'ccache'), CFLAGS='', CXXFLAGS='')
    if not args.resume:
        root.mkdir(parents=True, exist_ok=False)
    with destination_lock(root):
        build_locked(args, root, pixi, patches, expected, environment)


def build_locked(args, root, pixi, patches, expected, environment):
    source = root / 'source'
    record = root / 'build.json'
    if args.resume:
        report = json.loads(record.read_text())
        require_idle(report)
        if report.get('commit') != COMMIT or report.get('patches') != expected:
            raise ValueError('The retained build has different source/patch requirements.')
    else:
        report = {'status': 'preparing', 'commit': COMMIT, 'patches': expected,
                  'pixi': pixi, 'installer_sha256': digest(Path(__file__)), 'steps': []}

    report.setdefault('invocations', []).append({
        'installer_sha256': digest(Path(__file__)), 'resume': args.resume,
        'prepare_only': args.prepare_only, 'jobs': args.jobs,
    })

    def save():
        temporary = record.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(report, indent=2) + '\n')
        temporary.replace(record)

    def run(label, command):
        report['steps'].append({'label': label, 'command': list(map(str, command))})
        save()
        print(label, flush=True)
        with (root / 'build.log').open('ab') as log:
            log.write((json.dumps(list(map(str, command))) + '\n').encode()); log.flush()
            run_recorded(command, root, environment, log, report['steps'][-1], save)

    def inspect_source():
        head = subprocess.check_output(['git', '-C', source, 'rev-parse', 'HEAD'], text=True).strip()
        if head != COMMIT or digest(source / 'pixi.lock') != LOCK:
            raise ValueError('FreeCAD commit or locked SDK differs from qualification.')
        raw = subprocess.check_output(['git', '-C', source, 'diff', '--binary', 'HEAD'])
        return hashlib.sha256(raw).hexdigest()

    try:
        if not args.resume:
            run('Fetch pinned FreeCAD', ['git', 'clone', '--depth', '1', '--branch', '1.1.3',
                                        'https://github.com/FreeCAD/FreeCAD.git', source])
            inspect_source()
            for name in PATCHES:
                run('Check ' + name, ['git', '-C', source, 'apply', '--check', patches / name])
                run('Apply ' + name, ['git', '-C', source, 'apply', patches / name])
            report['source_diff_sha256'] = inspect_source()
            report['status'] = 'prepared'; save()
        elif inspect_source() != report.get('source_diff_sha256'):
            raise ValueError('Tracked FreeCAD sources changed since preparation.')
        report['status'] = 'prepared'
        report.pop('error', None)
        save()
        if args.prepare_only:
            print(f'Prepared sources only. Resume with --resume in {root}.')
            return
        report['status'] = 'building'; save()
        run('Fetch pinned submodules', ['git', '-C', source, 'submodule', 'update', '--init', '--recursive'])
        prefix = [pixi, 'run', '--locked', '--manifest-path', source / 'pixi.toml']
        run('Install locked SDK', [pixi, 'install', '--locked', '--manifest-path', source / 'pixi.toml'])
        run('Configure host', prefix + ['cmake', '--preset', 'conda-linux-release', '-S', source,
                                      '-B', source / 'build/release',
                                      f'-DCMAKE_JOB_POOLS=compile_jobs={args.jobs};link_jobs=1'])
        run('Build host', prefix + ['cmake', '--build', source / 'build/release', '--parallel', str(args.jobs)])
        if inspect_source() != report['source_diff_sha256']:
            raise ValueError('Tracked sources changed during the build.')
        binary = source / 'build/release/bin/FreeCAD'
        if not binary.is_file():
            raise ValueError('Build did not produce the FreeCAD executable.')
        launcher = root / 'FreeCAD-Vinkulum'
        lines = ['#!/bin/sh', 'set -eu', 'unset LD_PRELOAD PYTHONHOME PYTHONPATH LD_LIBRARY_PATH QT_PLUGIN_PATH QT_QPA_PLATFORM_PLUGIN_PATH']
        for key in ('PIXI_HOME', 'PIXI_CACHE_DIR'):
            lines.append(f'export {key}={shlex.quote(environment[key])}')
        lines.append('exec ' + shlex.join(list(map(str, prefix + [binary]))) + ' "$@"')
        launcher.write_text('\n'.join(lines) + '\n'); launcher.chmod(0o755)
        report.update(status='built', launcher=str(launcher), binary_sha256=digest(binary))
        save(); print(f'Built host: {launcher}\nNative GUI qualification must be run separately.')
    except Exception as error:
        report.update(status='failed', error=str(error)); save()
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--pixi', default='pixi')
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--prepare-only', action='store_true')
    options = parser.parse_args()
    try:
        build(options)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, str(error) + '\n')
