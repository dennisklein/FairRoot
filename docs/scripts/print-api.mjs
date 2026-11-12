// Prints the API pages and their member anchors, for links from the docs pages.
//   node scripts/print-api.mjs                     all types, namespaces and enums
//   node scripts/print-api.mjs FairRun             only names that contain "FairRun"
//   node scripts/print-api.mjs FairTask --members  with the member anchors
// Links in pages leave out the base: /api/FairTask/#Exec
import path from 'node:path';
import { getModel } from './api.mjs';

const docsDir = path.resolve(new URL('..', import.meta.url).pathname);
const m = getModel({ xmlDir: path.join(docsDir, 'build/doxygen/xml'), repoRoot: path.resolve(docsDir, '..'), base: '' });
const args = process.argv.slice(2);
const members = args.includes('--members');
const filter = args.find((a) => !a.startsWith('--'));
const show = (name) => !filter || name.includes(filter);

const out = [];
const list = (rec) => {
  const seen = new Set();
  for (const x of rec.visible ?? []) {
    if (!x.anchor || seen.has(x.anchor) || x.page) continue;
    seen.add(x.anchor);
    out.push(`    #${x.anchor}  ${x.prot !== 'public' ? `(${x.prot}) ` : ''}${x.sigText.replace(/\s+/g, ' ')}`);
  }
};
if (show('global')) {
  out.push(`global scope  ${m.global.url}`);
  if (members) list(m.global);
}
for (const ns of m.namespaces()) {
  if (!show(ns.name)) continue;
  out.push(`namespace ${ns.name}  ${ns.url}`);
  if (members) list(ns);
}
for (const t of m.types()) {
  if (!show(t.name)) continue;
  out.push(`${t.kind} ${t.name}  ${t.url}  (${t.file}:${t.line}, FairRoot::${t.lib?.library})`);
  if (members) list(t);
}
for (const { e } of m.namespaceEnums()) {
  if (!show(e.qualified)) continue;
  out.push(`enum ${e.qualified}  ${e.page}  (${e.file}:${e.line})`);
  if (members) for (const v of e.values) out.push(`    #${v.anchor}  ${v.name}${v.init ? ' ' + v.init : ''}`);
}
console.log(out.join('\n'));
