//package dihfind;

import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.LorentzVector;
import java.lang.Math;


public class EvParticle {

  public int PID;
  public double mass;
  public double p[];

  public static final int kX = 0;
  public static final int kY = 1;
  public static final int kZ = 2;


  protected LorentzVector lvec;



  public EvParticle(int pid) {
    this(pid,-1000,-1000,1000);
  }


  public EvParticle(int pid, double px, double py, double pz) {
    PID = pid;
    p = new double[3];
    p[kX] = px;
    p[kY] = py;
    p[kZ] = pz;
    mass = findMass(PID);
    lvec = new LorentzVector();
    setLvec();
  }

  public void reset() {
    p[kX] = -1000;
    p[kY] = -1000;
    p[kZ] = -1000;
    setLvec();
  }


  private void setLvec() { 
    lvec.setPxPyPzM(p[kX],p[kY],p[kZ],mass);
  }


  public double findMass(int pid) {
    for(PDG parti : PDG.values()) {
      if(PID==parti.pid()) return parti.mass();
    }
    System.out.println("WARNING: unrecognized PID%n");
    return -1000; // (if PID unrecognized)
  }
  

  public void setP(double px, double py, double pz) {
    p[kX] = px;
    p[kY] = py;
    p[kZ] = pz;
    setLvec();
  }

  private boolean momentumIsSet() {
    return p[kX]>-800 && p[kY]>-800 && p[kZ]>-800;
  }

  public int getPid() { return PID; }
  public LorentzVector getLvec() { return lvec; }
  public double getMass() { return mass; }
  public double getPx() { return p[kX]; }
  public double getPy() { return p[kY]; }
  public double getPz() { return p[kZ]; }

  public double getPt() {
    if(momentumIsSet())
      return Math.hypot(p[kX],p[kY]);
    else return -1000;
  }

  public double getMom() {
    if(momentumIsSet())
      return Math.sqrt( Math.pow(p[kX],2) + 
                        Math.pow(p[kY],2) + 
                        Math.pow(p[kZ],2) );
    else return -1000;
  }

  public double getE() { 
    if(momentumIsSet()) return lvec.e(); 
    else return -1000;
  }

  public double getPhi() {
    if(momentumIsSet()) return lvec.phi(); 
    else return -1000;
  }


  public void print() {
    System.out.format("p=(%.3f, %.3f, %.3f)  E=%.3f  m=%.5f%n",
      this.getPx(), this.getPy(), this.getPz(), this.getE(), this.getMass());
  };
      


}
