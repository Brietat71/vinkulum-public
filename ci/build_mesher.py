"""Build the pinned Linux Gmsh/OCCT 8 mesher in a user-owned directory.

Requires git, CMake, Ninja, a C++ compiler and the Linux/X11 development headers.
Existing source overrides must be clean checkouts of the exact pinned commits.
The resulting installation uses its absolute OCCT prefix; do not move it.
"""

import argparse
import hashlib
import json
import platform
import subprocess
from pathlib import Path

SOURCES = {
    "occt": (
        "https://github.com/Open-Cascade-SAS/OCCT.git",
        "b8f597c677811d1f9f4d8a97f5ae2825c0353a42",
    ),
    "gmsh": (
        "https://gitlab.onelab.info/gmsh/gmsh.git",
        "91b4154a2ba9865335a548b1c146a38d0dea7141",
    ),
}


def run(arguments, log):
    print(f"{log.name}: {arguments[0]}", flush=True)
    with log.open("a") as stream:
        stream.write(json.dumps(list(map(str, arguments))) + "\n")
        stream.flush()
        subprocess.run(
            list(map(str, arguments)),
            cwd=log.parent,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=True,
        )


def source(root, name, override):
    url, commit = SOURCES[name]
    folder = Path(override).resolve() if override else root / (name + "-source")
    log = root / (name + "-source.log")
    if not folder.exists():
        run(["git", "init", folder], log)
        run(["git", "-C", folder, "remote", "add", "origin", url], log)
        run(["git", "-C", folder, "fetch", "--depth", "1", "origin", commit], log)
        run(["git", "-C", folder, "checkout", "--detach", "FETCH_HEAD"], log)
    head = subprocess.check_output(
        ["git", "-C", str(folder), "rev-parse", "HEAD"], text=True
    ).strip()
    dirty = subprocess.check_output(
        ["git", "-C", str(folder), "status", "--porcelain"], text=True
    )
    if head != commit or dirty:
        raise ValueError(f"{folder} must be a clean checkout of {commit}.")
    return folder


def main(args):
    if platform.system() != "Linux":
        raise ValueError("This recipe currently qualifies Linux builds only.")
    if not 1 <= args.jobs <= 32:
        raise ValueError("Use 1–32 build jobs.")
    root = args.directory.resolve()
    root.mkdir(parents=True, exist_ok=True)
    occt, gmsh = (
        source(root, "occt", args.occt_source),
        source(root, "gmsh", args.gmsh_source),
    )
    occt_build, gmsh_build = root / "occt-build", root / "gmsh-build"
    occt_prefix, gmsh_prefix = root / "occt-install", root / "gmsh-install"
    run(
        [
            "cmake",
            "-S",
            occt,
            "-B",
            occt_build,
            "-G",
            "Ninja",
            "-DCMAKE_BUILD_TYPE=Release",
            f"-DINSTALL_DIR={occt_prefix}",
            "-DCMAKE_INSTALL_RPATH=$ORIGIN",
            "-DCMAKE_EXPORT_NO_PACKAGE_REGISTRY=ON",
            "-DBUILD_LIBRARY_TYPE=Shared",
            "-DBUILD_MODULE_Draw=OFF",
            "-DBUILD_MODULE_Visualization=OFF",
            "-DBUILD_MODULE_ApplicationFramework=ON",
            "-DBUILD_MODULE_DataExchange=ON",
            "-DUSE_FREETYPE=OFF",
            "-DUSE_TBB=OFF",
            "-DUSE_RAPIDJSON=OFF",
            "-DUSE_TCL=OFF",
        ],
        root / "occt-configure.log",
    )
    run(
        ["cmake", "--build", occt_build, "--parallel", str(args.jobs)],
        root / "occt-build.log",
    )
    run(["cmake", "--install", occt_build], root / "occt-install.log")
    run(
        [
            "cmake",
            "-S",
            gmsh,
            "-B",
            gmsh_build,
            "-G",
            "Ninja",
            "-DCMAKE_BUILD_TYPE=Release",
            f"-DCMAKE_INSTALL_PREFIX={gmsh_prefix}",
            f"-DCMAKE_PREFIX_PATH={occt_prefix}",
            f"-DOCC_INC={occt_prefix}/include/opencascade",
            "-DENABLE_OCC=ON",
            "-DENABLE_OCC_CAF=ON",
            "-DENABLE_FLTK=OFF",
            "-DENABLE_BUILD_LIB=OFF",
            "-DENABLE_BUILD_SHARED=ON",
            "-DENABLE_BUILD_DYNAMIC=ON",
            f"-DCMAKE_INSTALL_RPATH={occt_prefix}/lib",
        ],
        root / "gmsh-configure.log",
    )
    if (
        "Found OpenCASCADE version 8.0.1"
        not in (root / "gmsh-configure.log").read_text()
    ):
        raise ValueError("Gmsh did not configure against the pinned OCCT 8.0.1.")
    run(
        ["cmake", "--build", gmsh_build, "--parallel", str(args.jobs)],
        root / "gmsh-build.log",
    )
    run(["cmake", "--install", gmsh_build], root / "gmsh-install.log")
    executable = gmsh_prefix / "bin/gmsh"
    info = subprocess.check_output(
        [str(executable), "-info"], stderr=subprocess.STDOUT, text=True
    )
    if "OCC version   : 8.0.1" not in info:
        raise ValueError("The installed mesher does not report OCCT 8.0.1.")
    (root / "engine-info.txt").write_text(info)
    files = [
        executable,
        *sorted((occt_prefix / "lib").glob("*.so*")),
        *sorted((gmsh_prefix / "lib").glob("*.so*")),
    ]
    manifest = {
        "format": "vinkulum-mesher-build",
        "schema": 1,
        "sources": {
            key: {"url": url, "commit": commit}
            for key, (url, commit) in SOURCES.items()
        },
        "platform": platform.platform(),
        "jobs": args.jobs,
        "installed_files_sha256": {
            str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files
            if p.is_file()
        },
        "note": "Build artifact fingerprints, not an attestation of all libraries loaded at runtime. Not a relocatable desktop bundle.",
    }
    (root / "build-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Choose this executable in Studio: {executable}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--occt-source", type=Path)
    parser.add_argument("--gmsh-source", type=Path)
    main(parser.parse_args())
