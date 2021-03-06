NAME="llvm"
LLVM_VERSION="3.6.2"
PATCH_VERSION="2"

VERSION="${LLVM_VERSION}p${PATCH_VERSION}"
PACKAGE="$NAME-$VERSION"

FILE="$NAME-$LLVM_VERSION.src.tar.xz"
SHA256SUM="f60dc158bfda6822de167e87275848969f0558b3134892ff54fced87e4667b94"

FILE_CFE="cfe-$LLVM_VERSION.src.tar.xz"
SHA256SUM_CFE="ae9180466a23acb426d12444d866b266ff2289b266064d362462e44f8d4699f3"

export CFLAGS="-O3 -std=c99"
export CXXFLAGS="-O3 -std=c++0x"

pkg_download() {
	wget -nc "http://www.insieme-compiler.org/ext_libs/$FILE"
	echo "$SHA256SUM  $FILE" | sha256sum -c
	wget -nc "http://www.insieme-compiler.org/ext_libs/$FILE_CFE"
	echo "$SHA256SUM_CFE  $FILE_CFE" | sha256sum -c
}

pkg_extract() {
	tar xf "$FILE"
	mv "$NAME-$LLVM_VERSION.src" "$PACKAGE"
	tar xf "$FILE_CFE"
	mv "cfe-$LLVM_VERSION.src" "$PACKAGE/tools/clang"
}

pkg_configure() {
	./configure --prefix="$PREFIX/$PACKAGE" \
		--enable-assert=yes \
		--enable-debug-runtime=no \
		--enable-debug-symbols=no \
		--enable-optimized=yes \
		--enable-shared=yes \
		--enable-bindings=none
}

pkg_build() {
	make -j "$SLOTS" clang-only REQUIRES_RTTI=1
}

pkg_install() {
	make -j "$SLOTS" clang-only install
}

pkg_cleanup() {
	rm -rf "$PACKAGE" "$FILE" "$FILE_CFE"
}
