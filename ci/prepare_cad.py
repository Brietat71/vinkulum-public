"""Fetch hash-pinned upstream sources and apply the reviewed OCCT 8 patches.

Usage: python ci/prepare_cad.py /tmp/vinkulum-cad-sources
Then install both returned directories together with ./apps/studio[cad,test].
No dependency constraint is ignored and no third-party package is monkey-patched.
"""

import hashlib
import json
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

PACKAGES = (
    (
        "build123d",
        "0.11.1",
        "9647c3fc7e1ef11d9c335787a54cf26e9bbc2191186342fe9755ec24a480033b",
    ),
    (
        "ocpsvg",
        "0.6.0",
        "f08da4347cc90ecd3565395e9bda5746d46ab8aafd6a2681bb03a9c321b54039",
    ),
)


def main(destination):
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    patches = Path(__file__).parent / "patches"
    for name, version, checksum in PACKAGES:
        folder = destination / f"{name}-{version}"
        patch = patches / f"{name}-{version}-occt8.patch"
        marker = folder / ".vinkulum-patch-sha256"
        patch_hash = hashlib.sha256(patch.read_bytes()).hexdigest()
        if marker.exists() and marker.read_text() == patch_hash:
            print(folder)
            continue
        if folder.exists():
            raise RuntimeError(
                f"{folder} exists without the expected patch; choose an empty destination."
            )
        archive = destination / f"{name}-{version}.tar.gz"
        if not archive.exists():
            with urllib.request.urlopen(
                f"https://pypi.org/pypi/{name}/{version}/json", timeout=30
            ) as response:
                metadata = json.load(response)
            source = next(f for f in metadata["urls"] if f["packagetype"] == "sdist")
            if source["digests"]["sha256"] != checksum:
                raise RuntimeError(f"Unexpected upstream digest: {name}")
            urllib.request.urlretrieve(source["url"], archive)
        if hashlib.sha256(archive.read_bytes()).hexdigest() != checksum:
            raise RuntimeError(f"Corrupt source archive: {archive}")
        with tarfile.open(archive) as stream:
            stream.extractall(destination, filter="data")
        # Upstream mixes CRLF and LF; patch files use LF consistently.
        for path in folder.rglob("*"):
            if path.suffix in (".py", ".toml") and path.is_file():
                raw = path.read_bytes()
                if b"\r\n" in raw:
                    path.write_bytes(raw.replace(b"\r\n", b"\n"))
        subprocess.run(
            ["patch", "--batch", "-p1", "-i", str(patch.resolve())],
            cwd=folder,
            check=True,
        )
        marker.write_text(patch_hash)
        print(folder)


if __name__ == "__main__":
    main(sys.argv[1])
