import org.jlab.jnp.hipo4.io.*;
import org.jlab.jnp.hipo4.data.*;

import java.lang.Math;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class simpleAnalyze{

  protected static final int kE = 0;
  protected static final int kP = 1;
  protected static final int kM = 2;
  protected static final int N = 3;

  public static void main(String[] args) {

    if (args.length == 0) {
      System.out.println("ERROR: Please enter a hipo file as the first argument");
      System.exit(0);
    }
    
    HipoReader reader = new HipoReader();
    reader.open(args[0]);

    Event ev = new Event();
    //Bank particleBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
    Bank particleBank = new Bank(reader.getSchemaFactory().getSchema("MC::Particle"));
    Bank runconfigBank = new Bank(reader.getSchemaFactory().getSchema("RUN::config"));

    int nFound=0;
    int nTotal=0;


    int hadrons[] = {kP,kM};
    double En[] = new double[N];
    double Px[] = new double[N];
    double Py[] = new double[N];
    double Pz[] = new double[N];
    double EnTmp[] = new double[N];
    double PxTmp[] = new double[N];
    double PyTmp[] = new double[N];
    double PzTmp[] = new double[N];
    double Pt[] = new double[N];
    boolean found[] = new boolean[N];
    int pid,idx,evnum,npart;
    int h,i;

    double PartMass[] = new double[N];
    PartMass[kE] = 0.000511;
    PartMass[kP] = 0.139571;
    PartMass[kM] = 0.139571;

    int PartPid[] = new int[N];
    PartPid[kE] = 11;
    PartPid[kP] = 211;
    PartPid[kM] = -211;

    String outstr;

    String outfileName = "javaOut.dat";


    try {
      BufferedWriter outfile = new BufferedWriter(new FileWriter(outfileName));

      reader.getEvent(ev,0); // init
      while(reader.hasNext()==true) {
        
        for(h=0; h<N; h++) {
          En[h] = -1;
          found[h] = false;
        };

        reader.nextEvent(ev);

        ev.read(runconfigBank);
        evnum = runconfigBank.getInt("event",0);

        ev.read(particleBank);
        npart = particleBank.getRows();

        for(i=0; i<npart; i++) {

          pid = particleBank.getInt("pid",i);
          for(h=0; h<N; h++) {
            if(pid == PartPid[h]) {
              idx = h;

              PxTmp[idx] = particleBank.getFloat("px",i);
              PyTmp[idx] = particleBank.getFloat("py",i);
              PzTmp[idx] = particleBank.getFloat("pz",i);

              EnTmp[idx] = Math.sqrt( Math.pow(PxTmp[idx],2)
                                    + Math.pow(PyTmp[idx],2)
                                    + Math.pow(PzTmp[idx],2)
                                    + Math.pow(PartMass[idx],2)
                                    );

              if(EnTmp[idx] > En[idx]) {
                En[idx] = EnTmp[idx];
                Px[idx] = PxTmp[idx];
                Py[idx] = PyTmp[idx];
                Pz[idx] = PzTmp[idx];
                Pt[idx] = Math.sqrt( Math.pow(Px[idx],2) + Math.pow(Py[idx],2) );
                found[idx] = true;
              };

            };

          };
        };

        if( found[kE] && found[kP] && found[kM] ) {
          outstr = Integer.toString(evnum);
          for(int j : hadrons) {
            //outstr += String.format(" %.2f %.2f %.2f",Px[j],Py[j],Pz[j]);
            outstr += String.format(" %.2f %.2f",En[j],Pt[j]);
          };
          outstr += "\n";
          outfile.write(outstr);
          System.out.println(outstr);
          nFound++;
        };
        nTotal++;
      };

      outfile.close();
      System.out.println(outfileName + " created.");
    } catch(IOException e) { e.printStackTrace(); };
    System.out.println(nFound+" / "+nTotal+" events have e-,pi+,pi-\n");
  };
};
