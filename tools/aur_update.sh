#!/bin/bash
v=$(git tag -l --sort='-*taggerdate'|head -n 1)

read -n 1 -p "Updating to $v. Continue? [Y/N] " GOON
echo

if [ "$GOON" = 'Y' ] ; then
	wget "https://github.com/heitzmann/gdspy/archive/$v.tar.gz"
	sha256=$(sha256sum "$v.tar.gz" | cut -d' ' -f1)

	for i in '' 2; do
		git clone "ssh://aur@aur.archlinux.org/python${i}-gdspy.git"

		cd "python${i}-gdspy"
		sed -i -e "s|pkgver=.*$|pkgver=${v:1}|" \
			-e "s|pkgrel=.*$|pkgrel=1|" \
			-e "s|^.*sums=.*$|sha256sums=('$sha256')|" \
					 PKGBUILD
		makepkg --printsrcinfo > .SRCINFO

		git diff

		read -n 1 -p "Commit update? [Y/N] " GOON
		echo

		if [ "$GOON" = 'Y' ] ; then
			git commit -am "Update to $v"
			git push
		fi

		cd -
		rm -rf "python${i}-gdspy"
	done
	rm "$v.tar.gz"
fi
