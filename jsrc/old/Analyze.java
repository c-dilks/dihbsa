//package dihfind;

import NovelFitters.NovelBaseFitter;
import NovelFitters.MyParticle;
import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;

import java.util.ArrayList;
import java.util.List;

import org.jlab.clas.physics.*;
import org.jlab.groot.data.H2F;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.groot.fitter.ParallelSliceFitter;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.groot.data.GraphErrors;
import org.jlab.groot.data.TDirectory;
import org.jlab.groot.data.H1F;

import java.io.*;
import java.nio.file.*;
import java.nio.file.attribute.*;
import static java.nio.file.FileVisitResult.*;
import static java.nio.file.FileVisitOption.*;
import java.util.*;
import java.util.Arrays;
import java.lang.Math;

public class Analyze{
  
  static public boolean debug;
  static public boolean debugDIS;

  protected static H1F elePxDist;
  protected static H1F elePyDist;
  protected static H1F elePzDist;
  protected static H1F elePtDist;
  protected static H1F eleMomDist;
  protected static H2F elePxPy;

  protected static H1F numPions[];

  protected static H1F piPxDist[];
  protected static H1F piPyDist[];
  protected static H1F piPzDist[];
  protected static H1F piPtDist[];
  protected static H1F piMomDist[];
  protected static H2F piPxPy[];

  protected static H2F dihPtCorr;
  protected static H1F dihMassDist;
  protected static H1F deltaPhiDist;
  protected static H2F dihPhiCorr;
  protected static H1F dihRPhiDist;
  protected static H1F dihSPhiDist;

  protected static H1F Q2Dist;
  protected static H1F XDist;
  protected static H1F WDist;
  protected static H2F Q2vsX;
  protected static H2F Q2vsW;

  protected String ws;
  protected static boolean eleFound;
  protected static boolean pionFound[];
  protected boolean validDIS;
  protected boolean validDih;
  protected boolean validEvent;
  protected int nValid;
  protected int pionCount[];

  protected HipoDataSource reader;
  protected NovelBaseFitter novel_fitter;
  protected EventFilter filter;
  protected float torus;

  protected File outfile;
  protected FileWriter fw;
  protected BufferedWriter bw;

  protected static LorentzVector vQ; // q vector
  protected static LorentzVector vDihS; // sum of dihadron vectors
  protected static LorentzVector vDihR; // difference of dihadron vectors


  

  protected double eleMom,eleTmpMom;
  protected double pionMom,pionTmpMom;

 
  protected double dihMass;

  static protected double massE;
  static protected double massP;
  static protected double massPi;

  protected double Q2, X, Nu, W;

  protected MyParticle part;

  protected static EvParticle ele,eleTmp;
  protected static EvParticle beam;
  protected static EvParticle target;
  protected static EvParticle pion[];
  protected static EvParticle pionTmp[];

  protected static final int kP = 0;
  protected static final int kM = 1;
  protected static final String piName[] = {"pip","pim"};
  protected static final String piTitle[] = {"pi+","pi-"};
  private int sgn;

  protected static final double pi = 3.141593;


  EventData currentEvent;
  EventData currentMCEvent;

  public static void main(String[] args) {
    if (args.length == 0) {
      // exits program if input directory not specified
      System.out.println("ERROR: Please enter a hipo file as the first argument");
      System.exit(0);
    }

    Analyze an = new Analyze();
    an.debug = false;
    an.debugDIS = false;

    massE = PDG.Electron.mass();
    massP = PDG.Proton.mass();
    massPi = PDG.Pip.mass();

    // define beam and target
    // - force the beam energy at 10.6, since PhysicsEvent::beamParticle may have bug ???
    beam = new EvParticle(PDG.Electron.pid(), 0, 0, 10.6);
    target = new EvParticle(PDG.Proton.pid(), 0, 0, 0);
    /*
    beam = new EvParticle(PDG.Proton.pid(),
        physEvent.beamParticle().px(),
        physEvent.beamParticle().py(),
        physEvent.beamParticle().pz());
    */


    if(debug) {
      System.out.println("BEAM:");
      beam.getLvec().print();
      System.out.println("TARGET:");
      target.getLvec().print();
    }


    // define observables
    ele = new EvParticle(PDG.Electron.pid());
    eleTmp = new EvParticle(PDG.Electron.pid());
    pion = new EvParticle[2];
    pion[kP] = new EvParticle(PDG.Pip.pid());
    pion[kM] = new EvParticle(PDG.Pim.pid());
    pionTmp = new EvParticle[2];
    pionTmp[kP] = new EvParticle(PDG.Pip.pid());
    pionTmp[kM] = new EvParticle(PDG.Pim.pid());

    pionFound = new boolean[2];


    
    vQ = new LorentzVector();
    vDihS = new LorentzVector();
    vDihR = new LorentzVector();

    double ptb = 2.5;

    elePxDist = new H1F("elePxDist","electron px",100,-ptb,ptb);
    elePyDist = new H1F("elePyDist","electron py",100,-ptb,ptb);
    elePzDist = new H1F("elePzDist","electron pz",100,0,14);
    elePtDist = new H1F("elePtDist","electron pT",100,0,ptb);
    eleMomDist = new H1F("eleMomDist","electron Mom",100,0,14);
    elePxPy = new H2F("elePxPy","electron py vs px",100,-ptb,ptb,100,-ptb,ptb);

    numPions = new H1F[2];
    piPxDist = new H1F[2];
    piPyDist = new H1F[2];
    piPzDist = new H1F[2];
    piPtDist = new H1F[2];
    piMomDist = new H1F[2];
    piPxPy = new H2F[2];

    for(int pp=0; pp<2; pp++) {
      numPions[pp] = new H1F("num"+piName[pp]+"Dist",piTitle[pp]+" multiplicity",10,0,10);
      piPxDist[pp] = new H1F(piName[pp]+"PxDist",piTitle[pp]+" px",100,-ptb,ptb);
      piPyDist[pp] = new H1F(piName[pp]+"PyDist",piTitle[pp]+" py",100,-ptb,ptb);
      piPzDist[pp] = new H1F(piName[pp]+"PzDist",piTitle[pp]+" pz",100,0,14);
      piPtDist[pp] = new H1F(piName[pp]+"PtDist",piTitle[pp]+" pT",100,0,ptb);
      piMomDist[pp] = new H1F(piName[pp]+"MomDist",piTitle[pp]+" Mom",100,0,14);
      piPxPy[pp] = new H2F(piName[pp]+"PxPy",piTitle[pp]+" py vs px",100,-ptb,ptb,100,-ptb,ptb);
    }

    dihPtCorr = new H2F("dihPtCorr","hadron pT correlation",100,0,ptb,100,0,ptb);
    dihMassDist = new H1F("dihMassDist","dihadron invariant mass",200,0,4);
    deltaPhiDist = new H1F("deltaPhiDist","deltaPhi distribution",100,-pi,pi);
    dihRPhiDist = new H1F("dihRPhiDist","R phi distribution",100,-pi,pi);
    dihSPhiDist = new H1F("dihSPhiDist","S phi distribution",100,-pi,pi);
    dihPhiCorr = new H2F("dihPhiCorr","hadron phi correlation",100,-pi,pi,100,-pi,pi);


    float Q2min = 0;
    float Q2max = 12;
    float Xmin = 0;
    float Xmax = 1;
    float Wmin = 0;
    float Wmax = 5;
    int Nbins = 100;
    Q2Dist = new H1F("Q2Dist","Q^2 distribution",Nbins,Q2min,Q2max);
    XDist = new H1F("XDist","x distribution",Nbins,Xmin,Xmax);
    WDist = new H1F("WDist","W distribution",Nbins,Wmin,Wmax);
    Q2vsX = new H2F("Q2vsX","Q^2 vs. x",Nbins,Xmin,Xmax,Nbins,Q2min,Q2max);
    Q2vsW = new H2F("Q2vsW","Q^2 vs. W",Nbins,Wmin,Wmax,Nbins,Q2min,Q2max);



    an.doAnalysis(args);
    an.plot();
  }

  
  void doAnalysis(String[] args) {
    nValid = 0;

    reader = new HipoDataSource();

    novel_fitter = new NovelBaseFitter(10.6,false,false);

    filter = new EventFilter("11:+211:-211:X+:X-:Xn");

    Path singlefilename = Paths.get(args[0]);
    PathMatcher matcher = FileSystems.getDefault().getPathMatcher("glob:**.{hipo}");

    if (matcher.matches(singlefilename)) {
      readData("", args[0]);
    } else {
      File folder = new File(args[0]);
      File[] listOfFiles = folder.listFiles();
      for (int iF = 0; iF < listOfFiles.length; iF++) {
        if (listOfFiles[iF].isFile()) {
          System.out.println("File " + listOfFiles[iF].getName());
          Path filename = Paths.get(listOfFiles[iF].getName());
          if (matcher.matches(filename)) {
            readData(args[0], listOfFiles[iF].getName());
          }
        }
      }
    }

    System.out.format("%n%nnValid=%d%n%n",nValid);
  } // em doAnalysis


  void readData(String args0, String filename) {
    System.out.println("matched" + filename);
    reader.open(args0 + filename); // open hipo file

    pionCount = new int[2];

    // loop through events
    while (reader.hasEvent() == true) {

      if(debug || debugDIS) System.out.format("%n[+++] BEGIN EVENT [+++]%n");
      resetEventVars();

      HipoDataEvent hipoEvent = (HipoDataEvent) reader.getNextEvent();
      PhysicsEvent physEvent = novel_fitter.getPhysicsEvent(hipoEvent);


      HipoDataBank runConfig = (HipoDataBank) hipoEvent.getBank("RUN::config");
      torus=runConfig.getFloat("torus",0); // -1 is inbending, +1 is outbending ???
      // IT SEEMS THAT ALL FILES IN ~/fromHarut are inbending (torus==-1)!!!

      //if(torus>0) continue; // test inbending
      //if(torus<0) continue; // test outbending


      // event filter
      if(filter.isValid(physEvent) == true) {

        // loop over particles in event
        for (int i=0; i<physEvent.count(); i++) {
          part = (MyParticle) physEvent.getParticle(i);
          /*
             if(debug) {
             System.out.println("found pid "+ part.pid());
             System.out.println("p=("+part.px()+","+part.py()+","+part.pz()+")");
             }
             */


          // electron handler
          if(part.pid()==PDG.Electron.pid()) 
            analyzeElectron(part);
          if(part.pid()==PDG.Pip.pid() || 
             part.pid()==PDG.Pim.pid() ) 
            analyzePion(part);




        } // eo particle loop
      } // eo event filter



      validDih = analyzeDihadron();
      validDIS = computeDISkinematics();

      validEvent = eleFound && 
                   pionFound[kP] && 
                   pionFound[kM] && 
                   validDih &&
                   validDIS;


      if(validEvent) {
        fillPlots();
        nValid++;
      }


    } // eo event loop
  } // em readData


  /*
  // to do
  public void outputData(String str) {
    try {
      outfile = new File("outfile.dat");
      outfile.createNewFile();
      fw = new FileWriter(outfile.getAbsoluteFile());
      bw = new BufferedWriter(fw);
    } catch(IOException e) { e.printStackTrace(); }
    bw.close();
  }
  */

  
  
  public void fillPlots() {

    this.elePxDist.fill(ele.getPx());
    this.elePyDist.fill(ele.getPy());
    this.elePzDist.fill(ele.getPz());
    this.elePtDist.fill(ele.getPt());
    this.eleMomDist.fill(ele.getMom());
    this.elePxPy.fill(ele.getPx(),ele.getPy());
    if(debug) {
      System.out.println("FILL ELECTRON PLOTS WITH:");
      ele.print();
    }

    for(int pp=0; pp<2; pp++) {
      this.piPxDist[pp].fill(pion[pp].getPx());
      this.piPyDist[pp].fill(pion[pp].getPy());
      this.piPzDist[pp].fill(pion[pp].getPz());
      this.piPtDist[pp].fill(pion[pp].getPt());
      this.piMomDist[pp].fill(pion[pp].getMom());
      this.piPxPy[pp].fill(pion[pp].getPx(),pion[pp].getPy());
      this.numPions[pp].fill(pionCount[pp]);
    };


    this.dihPtCorr.fill(pion[kM].getPt(),pion[kP].getPt());

    double phi[] = new double[2];
    for(int pp=0; pp<2; pp++) {
      phi[pp] = correctAngle(pion[pp].getPhi());
    }
    double deltaPhi = correctAngle(phi[kP] - phi[kM]);
    double Rphi = correctAngle(vDihR.phi());
    double Sphi = correctAngle(vDihS.phi());

    this.dihPhiCorr.fill(pion[kM].getPhi(),pion[kP].getPhi());
    this.deltaPhiDist.fill(deltaPhi);
    this.dihRPhiDist.fill(Rphi);
    this.dihSPhiDist.fill(Sphi);
    this.dihMassDist.fill(dihMass);

    this.Q2Dist.fill(Q2);
    this.XDist.fill(X);
    this.Q2vsX.fill(X,Q2);
    this.Q2vsW.fill(W,Q2);
    this.WDist.fill(W);
    
    //outputData(String.format("%.3f %.3f",X,Q2)); // todo

    return;
  };

        
  public double correctAngle(double a) {
    while(a>pi) a -= 2*pi;
    while(a<-pi) a += 2*pi;
    return a;
  }

  public void analyzePion(MyParticle part) {
    if(part.pid() == PDG.Pip.pid()) sgn = kP;
    else if(part.pid() == PDG.Pim.pid()) sgn = kM;
    else {
      System.out.println("ERROR in Analyze.analyzePion: bad pion PID");
      return;
    }

    pionTmp[sgn].setP(part.px(), part.py(), part.pz());
    pionTmpMom = pionTmp[sgn].getMom();

    pionMom = pion[sgn].getMom();
    if(pionTmpMom > pionMom) {
      pion[sgn].setP(pionTmp[sgn].getPx(),pionTmp[sgn].getPy(),pionTmp[sgn].getPz());
    }

    if(debug) {
      System.out.format("PION sgn=%d",sgn);
      System.out.println("PION pionTmp:");
      pionTmp[sgn].print();
      System.out.println("PION pion:");
      pion[sgn].print();
    };

    pionCount[sgn]++;

    pionFound[sgn] = true;


    return;
  };


  public void analyzeElectron(MyParticle part) {
    eleTmp.setP(part.px(), part.py(), part.pz());
    eleTmpMom = eleTmp.getMom();

    //if(eleTmpMom<2.0) return; // momentum cut

    eleMom = ele.getMom();
    if(eleTmpMom > eleMom) {
      if(eleMom>0) System.out.println("found another electron with larger momentum");
      ele.setP(eleTmp.getPx(),eleTmp.getPy(),eleTmp.getPz());
    }

    if(debug) {
      System.out.println("ELECTRON eleTmp:");
      eleTmp.print();
      System.out.println("ELECTRON ele:");
      ele.print();
    };



    eleFound = true;
    return;
  }

  public boolean analyzeDihadron() {
    if(!pionFound[kP] || !pionFound[kM]) return false;
    
    vDihS.copy(pion[kP].getLvec());
    vDihS.add(pion[kM].getLvec());
    dihMass = vDihS.mass();

    vDihR.copy(pion[kP].getLvec());
    vDihR.sub(pion[kM].getLvec());




    return true;
  }


  public boolean computeDISkinematics() {
    if(!eleFound) return false;

    if(debugDIS) {
      System.out.println("--------------------");
      System.out.println("beam 4-momentum:");
      beam.getLvec().print();
      System.out.println("electron 4-momentum:");
      ele.getLvec().print();
    };

    vQ.copy(beam.getLvec());
    vQ.sub(ele.getLvec());

    Q2 = -1 * vQ.mass2();
    //Nu = target.getE() * vQ.e() / massP;
    Nu = beam.getE() - ele.getE();
    if(Nu<0) {
      System.out.println("WARNING: beam energy < electron energy");
      return false;
    }

    X = Q2 / ( 2 * massP * Nu );
    W = Math.sqrt(Math.pow(massP,2) + 2*massP*Nu - Q2);


    if(debugDIS) {
      System.out.println("q 4-momentum:");
      vQ.print();
      System.out.format("X=%.3f  Q2=%.3f  Nu=%.3f  W=%.3f%n",X,Q2,Nu,W);
      System.out.println("--------------------");
    };


    //if(W<2.0) return false; // W cut


    return true;
  };





  public void plot() {
    System.out.println("now plotting...");
    
    EmbeddedCanvas canv_ele = new EmbeddedCanvas();
    canv_ele.setSize(900,600);
    canv_ele.divide(3,2);
    canv_ele.cd(0);
    canv_ele.draw(elePxDist);
    canv_ele.cd(1);
    canv_ele.draw(elePyDist);
    canv_ele.cd(2);
    canv_ele.draw(elePzDist);
    canv_ele.cd(3);
    canv_ele.draw(elePxPy,"colz");
    canv_ele.cd(4);
    canv_ele.draw(elePtDist);
    canv_ele.cd(5);
    canv_ele.draw(eleMomDist);
    canv_ele.save("eleP.png");

    EmbeddedCanvas[] canv_pi = new EmbeddedCanvas[2];
    for(int pp=0; pp<2; pp++) {
      canv_pi[pp] = new EmbeddedCanvas();
      canv_pi[pp].setSize(900,600);
      canv_pi[pp].divide(4,2);
      canv_pi[pp].cd(0);
      canv_pi[pp].draw(piPxDist[pp]);
      canv_pi[pp].cd(1);
      canv_pi[pp].draw(piPyDist[pp]);
      canv_pi[pp].cd(2);
      canv_pi[pp].draw(piPzDist[pp]);
      canv_pi[pp].cd(3);
      canv_pi[pp].draw(piPxPy[pp],"colz");
      canv_pi[pp].cd(4);
      canv_pi[pp].draw(piPtDist[pp]);
      canv_pi[pp].cd(5);
      canv_pi[pp].draw(piMomDist[pp]);
      canv_pi[pp].cd(6);
      canv_pi[pp].draw(numPions[pp]);
      canv_pi[pp].save(piName[pp]+"P.png");
    }

    EmbeddedCanvas canv_corr = new EmbeddedCanvas();
    canv_corr.setSize(900,600);
    canv_corr.divide(3,2);
    canv_corr.cd(0);
    canv_corr.draw(dihPtCorr,"colz");
    canv_corr.cd(1);
    canv_corr.draw(dihMassDist);
    canv_corr.cd(2);
    canv_corr.draw(dihPhiCorr,"colz");
    canv_corr.cd(3);
    canv_corr.draw(deltaPhiDist);
    canv_corr.cd(4);
    canv_corr.draw(dihRPhiDist);
    canv_corr.cd(5);
    canv_corr.draw(dihSPhiDist);
    canv_corr.save("corr.png");

    EmbeddedCanvas canvDIS = new EmbeddedCanvas();
    canvDIS.setSize(900,600);
    canvDIS.divide(3,2);
    canvDIS.cd(0);
    canvDIS.draw(Q2Dist);
    canvDIS.cd(1);
    canvDIS.draw(XDist);
    canvDIS.cd(2);
    canvDIS.draw(WDist);
    canvDIS.cd(3);
    //canvDIS.getPad(3).getAxisX().setRange(0.001,1);
    //canvDIS.getPad(3).getAxisY().setRange(0.1,5);
    //canvDIS.getPad(3).getAxisX().setLog(true);
    //canvDIS.getPad(3).getAxisY().setLog(true);
    canvDIS.getPad(3).getAxisZ().setLog(true);
    canvDIS.draw(Q2vsX,"colz");
    canvDIS.cd(4);
    //canvDIS.getPad(4).getAxisX().setRange(0.001,1);
    //canvDIS.getPad(4).getAxisY().setRange(0.1,5);
    //canvDIS.getPad(4).getAxisX().setLog(true);
    //canvDIS.getPad(4).getAxisY().setLog(true);
    canvDIS.getPad(4).getAxisZ().setLog(true);
    canvDIS.draw(Q2vsW,"colz");
    canvDIS.save("dis.png");
  }




  public void resetEventVars() {
    eleFound = false;
    ele.reset();

    for(int pp=0; pp<2; pp++) {
      pionFound[pp] = false;
      pion[pp].reset();
      pionCount[pp] = 0;
    }

    dihMass=-1000;

    validDIS = false;
    validDih = false;
    validEvent = false;
  }

}


