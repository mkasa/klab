imgcat2
========

This tool displays images inline in your terminal. It is a modified version of
the `imgcat` script that ships with iTerm2
(https://iterm2.com/3.2/documentation-images.html), so it requires a terminal
that understands iTerm2's inline image protocol (iTerm2, WezTerm, etc.). It
also works through tmux/screen and over ssh.

The main difference from the original `imgcat` is that imgcat2 can control the
display size of the image. By default it scales images down to 35% of the
terminal width so that large images do not overflow the screen.

basic usage
------------
You can display one or more image files.

    imgcat2 foo.png

    imgcat2 a.png b.jpg c.gif

You can also pipe an image into it from stdin.

    cat foo.png | imgcat2

    curl -s https://example.com/foo.png | imgcat2

print
------
Give --print (or -p) to print the file name before the image.

    imgcat2 -p foo.png

width
------
Give --width (or -W) to set the display width. The default is 35%.

    imgcat2 -W 50% foo.png

The value follows the iTerm2 spec:

- N    : number of character cells (columns)
- Npx  : number of pixels (e.g. 300px)
- N%   : percentage of the terminal width (e.g. 50%)

height
-------
Give --height (or -H) to set the display height. By default the height is
chosen to keep the aspect ratio.

    imgcat2 -H 20 foo.png

The value follows the same spec as --width (cells, Npx, or N%).

You can combine both to fit an image into a box.

    imgcat2 -W 80% -H 40% foo.png

url
----
Give --url (or -u) to download an image from a URL and display it directly.
This requires curl.

    imgcat2 -u https://example.com/foo.png

help
-----
Give --help (or -h) to see the usage.

    imgcat2 --help
