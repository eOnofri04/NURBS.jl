# Markdown fast reference

_By Elia Onofri_

A markdown file is a `.md` or `.markdown` file.  
Mostly, Markdown is just regular text with a few non-alphabetic characters thrown in, like [\#] or [\*].

Moreover, like in Latex, in order to achive a line break you have even to let
two blank lines or a single space as last char of a line.

The following is a list of the most important and used commands:

## Headers
The headers are marked with some \# (1 to 6) at the beginning of a line.
The less \# you will put, the bigger would be the text.

## Text
In order to make *Italic text* you can mark a phrase with an `*` (or even a `_`)
at the beginning and another one at the end.

On the other hand, in order to make **Bold texts**, you have to use two `*` (or `_`).

You can even combine them in order to make \*\***Bold and \__Italic_\_ at
the same time**\*\*.

Last but not least, in order to have ~~Crossed Out Text~~, you have to use two
tilde symbols (`~`)

## Lists
You can make different type of lists by choosing a string to put at the beginning
of each rows they are made of. In order to make a second height list you can
separate then using spaces.

Kind of List | Beginning String
------------ | -----------------
Numbered lists | `1.`
Bullet list | `-` or `\*`
Marked List | `- [X]`
Unmarked List | `- [ ]`


Numbered List
1. First issue
1. Second Issue
  1. Second level first issue
  1. Second level second issue

Bullet List
- First level issue
- Second level issue
 - Second level issue
- First level issue

Marked List
- [X] Marked issue
- [ ] Unmarked issue

## Tables
In order to create a table you have to put:
 - on the first line all the Header divided by pipe symbols (`|`);
 - on the second line a line made of hyphens (`-`) divided by pype symbols;
 - on the seguent lines the entries divided by pipe symbols.
like this:

Header 1 \| Header 2  
\-\-\-\-\|\-\-\-\-  
1.1 \| 1.2  
2.1 \| 2.2  

in order to get this result:

Header 1 | Header 2
----|----
1.1 | 1.2
2.1 | 2.2

## Code
In order to highline some code you can use multiple solution:
- For inline code you can put it inside backsticks (\`);
- For new line code you can use four (4) spaces at the beginning of the line;
- For multiple line code you can use three backsticks in the a first row and three
  in the last one (the lines need to be blank). Morover you can specify a specific
  language for the code by putting it after the first three backsticks.

A sample of inline code: \``int i = 0;`\`.

A sample line code (four spaces):

    int i = 0;
   
A multiple line code:

\`\`\`
```
int i = 0;
while (i < n)
    System.out.println("Hello World")
```
\`\`\`

A specific language multiple line code:

\`\`\`java
```java
int i = 0;
while (i < n)
    System.out.println("Hello World")
```
\`\`\`

## URLs
The markdown language automatically identify URLs like the following one:
http://gitlhub.com

However you can associate a URL with a String. In order to do so you have to put
the name you want to between box brackets and the URL between round brackets:

For esample you can write `[GitHub](http://github.com)` in order to obtain the
following link form: [GitHub](http://github.com).

## Extras

If you want to make a quotation you can use the `>` at the beginning of a row,
like this:

> \> In the whole of History
> \> Never before
> \> was so much owed
> \> by so many to so few
>
> \> W. Churchil

