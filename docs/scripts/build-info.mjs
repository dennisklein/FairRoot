// Version and commit of the FairRoot checkout that the site is built from.
import { execFileSync } from 'node:child_process';
import { readFileSync } from 'node:fs';
import path from 'node:path';
import { repoUrl } from './site.mjs';

function git(cwd, ...args) {
  try {
    return execFileSync('git', args, { cwd, encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] }).trim();
  } catch {
    return '';
  }
}

/** The version in `project(FairRoot VERSION x.y.z)` of the top CMakeLists.txt. */
export function projectVersion(repoRoot) {
  const text = readFileSync(path.join(repoRoot, 'CMakeLists.txt'), 'utf8');
  const m = /project\(\s*FairRoot\s+VERSION\s+(\d+)\.(\d+)\.(\d+)/i.exec(text);
  return m ? { major: Number(m[1]), minor: Number(m[2]), patch: Number(m[3]) } : null;
}

export function buildInfo(cwd = process.cwd()) {
  const commit = process.env.GITHUB_SHA || git(cwd, 'rev-parse', 'HEAD') || 'unknown';
  const commitShort = commit.slice(0, 7);
  const date = git(cwd, 'show', '-s', '--format=%cs', 'HEAD');
  const branch = process.env.GITHUB_HEAD_REF || process.env.GITHUB_REF_NAME || git(cwd, 'rev-parse', '--abbrev-ref', 'HEAD');
  const v = projectVersion(cwd);
  const version = v ? `${v.major}.${v.minor}.${v.patch}` : '';
  // A release build has the tag of its version on HEAD. Other builds show the
  // release series of the CMake version, with a "dev" mark.
  const tags = git(cwd, 'tag', '--points-at', 'HEAD').split('\n').filter(Boolean);
  const release = version && tags.includes(`v${version}`) ? `v${version}` : '';
  const label = release || (v ? `v${v.major}.${v.minor}` : 'dev');
  return {
    commit,
    commitShort,
    date,
    branch,
    version,
    release,
    label,
    unreleased: !release,
    repo: repoUrl,
  };
}
