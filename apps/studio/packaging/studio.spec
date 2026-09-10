from pathlib import Path
from importlib.metadata import version
import os
import sys
import tomllib
from packaging.version import Version
from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs, collect_submodules, copy_metadata
import OCP

assert int(OCP.__version__.split('.')[0]) >= 8, 'OCCT 8 minimum required for Studio CAD'

root = Path(SPECPATH).parents[2]
studio_version = tomllib.loads((root / 'apps/studio/pyproject.toml').read_text())['project']['version']
assert version('vinkulum-studio') == studio_version, 'Reinstall Studio before packaging'
datas = collect_data_files('vinkulum')
if sys.platform == 'darwin':
    payload = Path(os.environ['VINKULUM_MACOS_PAYLOAD'])
    datas += [(str(payload / 'Examples'), 'Examples'), (str(payload / 'build-info.json'), '.')]
datas += [(str(Path(SPECPATH) / 'licenses'), 'licenses')]
datas += [(str(root / 'ci' / 'patches'), 'licenses/cad-patches')]
datas += copy_metadata('build123d', recursive=True)
for package in ('OCP', 'build123d', 'ocpsvg', 'ocp_gordon', 'lib3mf'):
    datas += collect_data_files(package)
for distribution in ('vinkulum', 'vinkulum-studio', 'PySide6', 'PySide6_Essentials', 'shiboken6', 'vtk', 'numpy', 'pyinstaller', 'cadquery-ocp-novtk', 'cadquery-ocp-proxy', 'build123d', 'ocpsvg', 'ocp_gordon', 'lib3mf'):
    datas += copy_metadata(distribution)
a = Analysis([str(Path(SPECPATH) / 'launcher.py')],
    pathex=[str(root / 'apps/studio')], datas=datas, binaries=collect_dynamic_libs('lib3mf'),
    hiddenimports=['vinkulum._vinkulum', 'vtkmodules.vtkRenderingOpenGL2', 'vtkmodules.vtkInteractionStyle', 'vinkulum_studio.cad_worker'] + collect_submodules('OCP'),
    excludes=['PyQt5', 'PyQt6', 'PySide2', 'tkinter'], noarchive=False)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name='Vinkulum Studio',
          console=False, target_arch='arm64' if sys.platform == 'darwin' else None,
          codesign_identity=None)
coll = COLLECT(exe, a.binaries, a.datas, name='Vinkulum Studio')
if sys.platform == 'darwin':
    app = BUNDLE(coll, name='Vinkulum Studio.app', bundle_identifier='org.vinkulum.studio',
                 version=Version(studio_version).base_version,
                 info_plist={'NSHighResolutionCapable': True,
                 'LSMinimumSystemVersion': '14.0',
                 'CFBundleShortVersionString': Version(studio_version).base_version,
                 'VinkulumStudioVersion': studio_version})
