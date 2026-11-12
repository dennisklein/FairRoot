// Groups the members of a class by their name, for the filter on the API pages.

export const FAMILY_ORDER = [
  'Construction',
  'Types',
  'Lifecycle',
  'Get',
  'Set',
  'Register',
  'Create',
  'Run',
  'Other',
  'Operators',
  'Constants',
  'Data',
  'Friends',
];

// The hooks that FairRoot calls on tasks, modules, detectors, generators and runs.
const LIFECYCLE =
  /^(Init|InitTask|ReInit|ReInitTask|SetParContainers|SetParTask|InitParContainers|Exec|ExecuteTask|ExecMQ|InitMQ|Finish|FinishTask|FinishEvent|FinishRun|FinishPrimary|BeginEvent|BeginPrimary|BeginWorkerRun|FinishWorkerRun|PreTrack|PostTrack|EndOfEvent|Initialize|ProcessHits|ConstructGeometry|ConstructOpGeometry|ConstructRootGeometry|ConstructASCIIGeometry|ConstructGDMLGeometry|ReadEvent|Reset|Stepping|PreInit|PostInit|Setup|SetupPostInit|Close|Open|Fill)$/;

export function family(m, owner) {
  const n = m.name;
  if (m.kind === 'friend') return 'Friends';
  if (m.kind === 'typedef' || m.kind === 'enum') return 'Types';
  if (m.kind === 'variable') return m.static ? 'Constants' : 'Data';
  if (m.kind === 'function') {
    if (n === owner.short || n === `~${owner.short}` || /^operator\s*=$/.test(n) || n === 'Instance' || /^Clone/.test(n)) return 'Construction';
    if (/^operator/.test(n)) return 'Operators';
    if (LIFECYCLE.test(n)) return 'Lifecycle';
    if (/^(Get|Is|Has|Are|Can|Was|Use|Check)[A-Z_]/.test(n) || /^(get|is|has)[A-Z_]/.test(n)) return 'Get';
    if (/^(Set|Update|Change|Reset|Add|Remove|Activate|Enable|Disable)[A-Z_]/.test(n) || /^(set|add)[A-Z_]/.test(n)) return 'Set';
    if (/^(Register|Unregister)/.test(n)) return 'Register';
    if (/^(New|Create|Make|Build)[A-Z_]/.test(n) || /^(create|make)[A-Z_]/.test(n)) return 'Create';
    if (/^(Run|Process|Read|Write|Store|Print|Draw|Terminate)/.test(n)) return 'Run';
  }
  return 'Other';
}
