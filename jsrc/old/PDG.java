public enum PDG {
  Pip(211, 0.13957),
  Pim(-211, 0.13957),
  Kap(321, 0.493677),
  Kam(321, 0.493677),
  Proton(2212, 0.938272),
  Electron(11, 0.000511);
  private final int mPid;
  private final double mMass;

  PDG(int pid, double mass) {
    this.mPid = pid;
    this.mMass = mass;
  }
  public int pid() { return mPid; }
  public double mass() { return mMass; }
}
