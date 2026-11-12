// Shared access to the API model for the generated pages.
import { readFileSync } from 'node:fs';
import path from 'node:path';
import { paths, build } from 'virtual:fairroot/data';
import { getModel, LIBRARY_ORDER, SOURCE_DIRS } from '../../scripts/api.mjs';

export const base = import.meta.env.BASE_URL.replace(/\/$/, '');

export function model() {
  return getModel({ xmlDir: paths.xmlDir, repoRoot: paths.repoRoot, base });
}

let theme;
export function codeTheme() {
  theme ??= JSON.parse(readFileSync(path.join(paths.docsDir, 'src/styles/code-theme.json'), 'utf8'));
  return theme;
}

/** A path without the site base, for Starlight sidebar links (Starlight adds the base). */
export const noBase = (href) => (href.startsWith(base + '/') ? href.slice(base.length) : href);

/** GitHub link to a file and line at the commit of the build. */
export function githubUrl(file, line) {
  const ref = build.commit && build.commit !== 'unknown' ? build.commit : 'dev';
  return `${build.repo}/blob/${ref}/${file}${line ? `#L${line}` : ''}`;
}

const KIND_LETTER = { class: 'C', struct: 'S', union: 'U', enum: 'E', namespace: 'N' };
const byName = (a, b) => a.localeCompare(b, 'en', { sensitivity: 'base' });
const libRank = (lib) => {
  const i = LIBRARY_ORDER.indexOf(lib);
  return i < 0 ? LIBRARY_ORDER.length : i;
};

/**
 * The public types and enums, grouped by library (the CMake target FairRoot::<Library>)
 * and by the directory in the library. Each group: { library, sub, id, entries }.
 */
export function apiGroups() {
  const m = model();
  const groups = new Map();
  const add = (lib, entry) => {
    const library = lib?.library ?? 'Other';
    const sub = lib?.sub ?? '';
    const key = `${library}/${sub}`;
    if (!groups.has(key)) {
      groups.set(key, { library, sub, module: lib?.module ?? '', id: `${library}${sub ? '-' + sub : ''}`.toLowerCase().replace(/[^a-z0-9]+/g, '-'), entries: [] });
    }
    groups.get(key).entries.push(entry);
  };
  for (const t of m.types()) {
    // Nested types follow their parent on the parent page; the sidebar shows them indented.
    add(t.lib, { name: t.name, short: t.short, url: t.url, kind: t.kind, brief: t.brief, nested: t.scope && m.compoundByName(t.scope)?.kind !== 'namespace' });
  }
  for (const { e } of m.namespaceEnums()) {
    add(m.libraryOf(e.file), { name: e.qualified, short: e.name, url: e.page, kind: 'enum', strong: e.strong, brief: e.brief });
  }
  const sortKey = (e) => e.name.split('::').map((s) => s.toLowerCase());
  const list = [...groups.values()];
  for (const g of list) {
    g.entries.sort((a, b) => {
      const ka = sortKey(a);
      const kb = sortKey(b);
      for (let i = 0; i < Math.min(ka.length, kb.length); i++) if (ka[i] !== kb[i]) return ka[i] < kb[i] ? -1 : 1;
      return ka.length - kb.length;
    });
  }
  return list.sort((a, b) => libRank(a.library) - libRank(b.library) || byName(a.library, b.library) || byName(a.sub, b.sub));
}

/** Sidebar for the API pages: one group per library, with a sub-group per directory. */
export function apiSidebar() {
  const m = model();
  const libs = new Map();
  for (const g of apiGroups()) {
    if (!libs.has(g.library)) libs.set(g.library, []);
    libs.get(g.library).push(g);
  }
  const item = (e) => ({
    label: e.short,
    link: noBase(e.url),
    attrs: { 'data-kind': KIND_LETTER[e.kind] ?? 'C', class: e.nested ? 'fr-mono fr-nested' : 'fr-mono', title: e.name },
  });
  const out = [
    { label: 'API reference', link: '/api/' },
    { label: 'Global scope', link: noBase(m.global.url), attrs: { 'data-kind': 'N' } },
    ...m.namespaces()
      .sort((a, b) => byName(a.name, b.name))
      .map((ns) => ({ label: ns.name, link: noBase(ns.url), attrs: { 'data-kind': 'N', class: 'fr-mono' } })),
  ];
  for (const [library, groups] of libs) {
    const items = [];
    for (const g of groups) {
      if (!g.sub) items.push(...g.entries.map(item));
      else items.push({ label: `${g.sub}/`, collapsed: true, items: g.entries.map(item) });
    }
    // A library with one directory needs no sub-group.
    const flat = items.length === 1 && items[0].items ? items[0].items : items;
    out.push({ label: library, collapsed: library !== 'Base', items: flat });
  }
  return out;
}

/** The id of a directory on the source overview, for example "fairroot-base-steer". */
export const dirId = (dir) => dir.replace(/[^A-Za-z0-9]+/g, '-');

/**
 * Sidebar for the source pages. It lists the directories of each top folder, and the
 * files of the directory of the current page. A directory links to its row on the
 * source overview, so that the sidebar stays short.
 */
export function sourceSidebar(current) {
  const m = model();
  const dirs = new Map();
  for (const f of m.sourceFiles) {
    const dir = path.posix.dirname(f);
    if (!dirs.has(dir)) dirs.set(dir, []);
    dirs.get(dir).push(f);
  }
  const currentDir = current ? path.posix.dirname(current) : null;
  const out = [{ label: 'Source browser', link: '/source/' }];
  for (const top of SOURCE_DIRS) {
    const items = [];
    for (const [dir, list] of dirs) {
      if (dir !== top && !dir.startsWith(top + '/')) continue;
      const label = dir === top ? `${top}/` : `${dir.slice(top.length + 1)}/`;
      if (dir === currentDir) {
        items.push({
          label,
          collapsed: false,
          items: list.map((f) => ({ label: path.posix.basename(f), link: `/source/${f}/`, attrs: { class: 'fr-mono' } })),
        });
      } else {
        items.push({ label, link: `/source/#${dirId(dir)}`, attrs: { class: 'fr-mono fr-dir' } });
      }
    }
    out.push({ label: `${top}/`, collapsed: !(currentDir && (currentDir === top || currentDir.startsWith(top + '/'))), items });
  }
  return out;
}

/** Docs pages that name a file in their `code` front matter. A page that names Foo.cxx also
 *  describes Foo.h, and the other way round. */
export function guidesForFile(entries, file) {
  const stem = (f) => f.replace(/\.(h|cxx)$/, '');
  return entries
    .filter((e) => (e.data.code ?? []).some((c) => c === file || stem(c) === stem(file)))
    .map((e) => ({ label: e.data.title, link: `/${e.id}/` }));
}
