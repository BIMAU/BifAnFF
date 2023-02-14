---
layout: page
title: Instructions
permalink: /instructions/
---

This page is for the maintainers of these pages.

In general, to generate web pages one needs to define the style and define content. In github, jekyll is used to generate web pages and one can also define the style in it.
For the content one uses a markdown language. We chose in "_config.yml" to use kramdown. The idea of a markdown language is that the ascii file also looks already structured and is
well readable.


We created these pages by creating on a local desktop a new branch of BifAnFF called docs and these were filled by a [creation call to jekyll](https://docs.github.com/en/pages/setting-up-a-github-pages-site-with-jekyll/creating-a-github-pages-site-with-jekyll).
This generates example files which can be pushed to github and will display the contents, after some delay (15 minutes), on <https://bimau.github.io/BifAnFF>.
These files can be amended and also be viewed on the local desktop when a jekyll server has been started. 

Each page consists of [style instructions for jekyll](https://jekyllrb.com/docs/) at the beginning between two dashed lines "------------". Next the content starts which is written in [kramdown](https://kramdown.gettalong.org/syntax.html).

Frontmatter content appearing on each page is defined in "_config.yml". The content for the main page can be found "index.md" and one may define more ".md" files by copying "about.md". These appear as links on all pages.
Next there is a directory "_posts" containing posts which by copying, and adjusting the name, can be turned into another post.
