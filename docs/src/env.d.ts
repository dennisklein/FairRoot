declare module 'virtual:fairroot/data' {
  /** Version and commit of the FairRoot checkout that the site is built from. */
  export const build: {
    commit: string;
    commitShort: string;
    date: string;
    branch: string;
    version: string;
    release: string;
    label: string;
    unreleased: boolean;
    repo: string;
  };
  /** Absolute paths, resolved when the site build starts. */
  export const paths: {
    docsDir: string;
    repoRoot: string;
    xmlDir: string;
  };
}
