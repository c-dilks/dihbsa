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
    Bank particleBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
    Bank runconfigBank = new Bank(reader.getSchemaFactory().getSchema("RUN::config"));

    double En[] = new double[N];
    double EnTmp[] = new double[N];
    double Pt[] = new double[N];
    double Px[] = new double[N];
    double Py[] = new double[N];
    double Pz[] = new double[N];
    double EnMax[] = new double[N];
    boolean found[] = new boolean[N];
    int npart;
    int pid;
    int idx;
    int evnum;

    double Mass[] = new double[N];
    Mass[kE] = 0.000511;
    Mass[kP] = 0.139571;
    Mass[kM] = 0.139571;

    String outstr;

    String outfileName = "javaOut.dat";


    try {
      BufferedWriter outfile = new BufferedWriter(new FileWriter(outfileName));

      reader.getEvent(ev,0); // init
      while(reader.hasNext()==true) {

        reader.nextEvent(ev);

        ev.read(runconfigBank);
        evnum = runconfigBank.getInt("event",0);

        ev.read(particleBank);
        npart = particleBank.getRows();

        if(npart>0) {

          for(int h=0; h<N; h++) {
            En[h] = -1;
            Pt[h] = -1;
            found[h] = false;
          };
          idx = -1;

          for(int i=0; i<npart; i++) {

            pid = particleBank.getInt("pid",i);
            if(pid == 11) idx = kE;
            else if(pid==211) idx = kP;
            else if(pid==-211) idx = kM;

            if(idx>=0) {
              Px[idx] = particleBank.getFloat("px",i);
              Py[idx] = particleBank.getFloat("py",i);
              Pz[idx] = particleBank.getFloat("pz",i);

              EnTmp[idx] = Math.sqrt(
                Math.pow(Px[idx],2) + Math.pow(Py[idx],2) + Math.pow(Pz[idx],2)
                - Math.pow(Mass[idx],2));

              if(EnTmp[idx] > En[idx]) {
                En[idx] = EnTmp[idx];
                Pt[idx] = Math.sqrt( Math.pow(Px[idx],2) + Math.pow(Py[idx],2) );
                found[idx] = true;
              };

            };
          };

          if( found[kE] && found[kP] && found[kM] ) {
            outstr = Integer.toString(evnum);
            for(int h=0; h<N; h++) 
              outstr += String.format(" %.2f %.2f",En[h],Pt[h]);
            outstr += "\n";
            outfile.write(outstr);
            //System.out.println(outstr);
          };
        };
      };

      outfile.close();
      System.out.println(outfileName + " created.");
    } catch(IOException e) { e.printStackTrace(); };
  };
};
