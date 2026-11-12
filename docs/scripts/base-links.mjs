// Sätteri hast plugin: prefixes root-relative links and images in Markdown and MDX
// with the site base, so that pages can link to `/concepts/tasks/`
// instead of `/FairRoot/concepts/tasks/`.
// Links to `/repo/<path>` go to the file or folder in the repository on GitHub:
// `/repo/CMakeLists.txt` becomes `https://github.com/<owner>/FairRoot/blob/dev/CMakeLists.txt`.
// A line anchor works too: `/repo/CHANGELOG.md#L10-L20` (Markdown files show as plain text then).
import { repoUrl, branch } from './site.mjs';

function repoLink(url) {
  const [p, hash] = url.slice('/repo/'.length).split('#');
  const plain = hash && /^L\d/.test(hash) && /\.(md|markdown)$/i.test(p) ? '?plain=1' : '';
  return `${repoUrl}/blob/${branch}/${p}${plain}${hash ? `#${hash}` : ''}`;
}

function prefix(url, base) {
  if (typeof url !== 'string' || !url.startsWith('/') || url.startsWith('//')) return url;
  if (url.startsWith('/repo/')) return repoLink(url);
  if (url === base || url.startsWith(base + '/')) return url;
  return base + url;
}

export default function baseLinks({ base }) {
  const b = base.replace(/\/$/, '');
  return {
    name: 'fairroot-base-links',
    element: {
      filter: ['a', 'img'],
      visit(node, ctx) {
        if (!b) return;
        const key = node.tagName === 'a' ? 'href' : 'src';
        const value = node.properties?.[key];
        const next = prefix(value, b);
        if (next !== value) ctx.setProperty(node, key, next);
      },
    },
  };
}
