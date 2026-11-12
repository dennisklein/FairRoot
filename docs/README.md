# FairRoot documentation site

This folder is the source of the FairRoot documentation site on GitHub Pages. The
site uses [Astro](https://astro.build/) with [Starlight](https://starlight.astro.build/).

## Requirements

* Node.js 22.23.4 or newer
* Doxygen 1.9 or newer. The API reference and the source browser use its XML output.
* git. The site shows the commit and the branch.

## Build

```sh
cd docs
npm ci
npm run dev       # local server at http://localhost:4321/FairRoot/, it reloads after changes
npm run build     # writes the site to dist/
npm run check     # checks all links in dist/
npm run preview   # serves dist/
```

The build runs Doxygen on `fairroot/`, `examples/` and `templates/`, and writes the
XML to `build/doxygen/`. The next builds run Doxygen again only when a C++ file in
these folders is newer than the XML.

The site URL comes from `GITHUB_REPOSITORY` (see `scripts/site.mjs`). A fork gets its
own site with links to its own repository.

## Contents

| Path | Contents |
| --- | --- |
| `src/content/docs/` | The pages, in MDX. Each folder is a group of the sidebar. |
| `src/pages/api/` | The API reference, from the Doxygen XML. The pages are grouped by library. |
| `src/pages/source/` | The source browser, from the Doxygen XML |
| `src/pages/404.astro` | The 404 page. It sends the URLs of the old Doxygen site to the new pages. |
| `src/components/` | Components for the pages: `Snippet` and the diagrams |
| `scripts/` | The data for the generated pages, and the checks |

## Write a page

* Each page starts with frontmatter: `title`, `description` and `sidebar.order`.
  `code` lists the source files that the page describes. The source browser links
  back from these files to the page.
* Show code of the repository with `<Snippet>`, not with a copy. Then the samples
  are code that CI compiles, and the build fails when the code changes and the
  sample does not agree. `src/components/Snippet.astro` shows the usage.
* Write links root-relative and without the base: `/concepts/tasks/`,
  `/api/FairTask/`, `/source/fairroot/base/steer/FairTask.h/#L36`.
* Run `node scripts/check-pages.mjs` before you build. It checks the frontmatter, the
  MDX syntax, the snippets and all links in a few seconds.
  `node scripts/print-api.mjs FairTask --members` prints the API pages and anchors
  of a class.
* Write in the style of ASD-STE100 (Simplified Technical English): short sentences,
  one instruction for each sentence, the active voice.
