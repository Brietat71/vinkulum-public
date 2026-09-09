fn main() {
    for file in [
        "abi.cpp",
        "inertie_binaire_native.cpp",
        "comparaison_masse_native.cpp",
    ] {
        println!("cargo:rerun-if-changed=src/certificat/{file}");
    }
    let mut build = cc::Build::new();
    build
        .cpp(true)
        .std("c++17")
        .opt_level(3)
        .file("src/certificat/abi.cpp");
    let compiler = build.get_compiler();
    if compiler.is_like_msvc() {
        build.flag("/fp:strict");
    } else if compiler.is_like_gnu() || compiler.is_like_clang() {
        build.flags(["-fno-fast-math", "-ffp-contract=off", "-frounding-math"]);
    } else {
        panic!("Les certificats requièrent GCC, Clang ou MSVC avec arithmétique stricte");
    }
    let compiler = build.get_compiler();
    println!(
        "cargo:rustc-env=VINKULUM_CERTIFICAT_COMPILATEUR={}",
        compiler.path().display()
    );
    println!(
        "cargo:rustc-env=VINKULUM_CERTIFICAT_OPTIONS={:?}",
        compiler.args()
    );
    build.compile("vinkulum_certificat");
}
