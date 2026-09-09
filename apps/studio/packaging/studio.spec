from pathlib import Path
from PyInstaller.utils.hooks import collect_data_files, copy_metadata

root = Path(SPECPATH).parents[2]
datas = collect_data_files('vinkulum')
datas += [(str(Path(SPECPATH) / 'licenses'), 'licenses')]
for distribution in ('vinkulum', 'vinkulum-studio', 'PySide6', 'PySide6_Essentials', 'shiboken6', 'vtk', 'numpy', 'pyinstaller'):
    datas += copy_metadata(distribution)
a = Analysis([str(Path(SPECPATH) / 'launcher.py')],
    pathex=[str(root / 'apps/studio')], datas=datas,
    hiddenimports=['vinkulum._vinkulum', 'vtkmodules.vtkRenderingOpenGL2', 'vtkmodules.vtkInteractionStyle'],
    excludes=['PyQt5', 'PyQt6', 'PySide2', 'tkinter'], noarchive=False)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name='Vinkulum Studio',
          console=False, target_arch='arm64', codesign_identity=None)
coll = COLLECT(exe, a.binaries, a.datas, name='Vinkulum Studio')
app = BUNDLE(coll, name='Vinkulum Studio.app', bundle_identifier='org.vinkulum.studio',
             version='0.3.0', info_plist={'NSHighResolutionCapable': True,
             'LSMinimumSystemVersion': '14.0', 'CFBundleShortVersionString': '0.3.0'})
