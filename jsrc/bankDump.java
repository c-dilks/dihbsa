import org.jlab.jnp.hipo4.io.*;
import org.jlab.jnp.hipo4.data.*;

import java.lang.Math;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class bankDump{
  
  static int evnum;

  public static void main(String[] args) {

    if (args.length == 0) {
      System.out.println("ERROR: Please enter a hipo file as the first argument");
      System.exit(0);
    }
    
    HipoReader reader = new HipoReader();
    reader.open(args[0]);

    int whichEv = 0;
    if(args.length == 2) whichEv = Integer.parseInt(args[1]);


    Event ev = new Event();
    Bank runconfigBank = new Bank(reader.getSchemaFactory().getSchema("RUN::config"));
    Bank parBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
    //Bank genBank = new Bank(reader.getSchemaFactory().getSchema("MC::Particle"));
    //Bank lundBank = new Bank(reader.getSchemaFactory().getSchema("MC::Lund"));
    Bank calBank = new Bank(reader.getSchemaFactory().getSchema("REC::Calorimeter"));
    Bank trkBank = new Bank(reader.getSchemaFactory().getSchema("REC::Track"));
    Bank trajBank = new Bank(reader.getSchemaFactory().getSchema("REC::Traj"));

    reader.getEvent(ev,0); // init

    int npart;

    while(reader.hasNext()==true) {

      reader.nextEvent(ev);

      ev.read(runconfigBank);
      evnum = runconfigBank.getInt("event",0);

      if(whichEv==0 || whichEv==evnum) {
        System.out.println("\n-----");
        System.out.println("EVNUM = "+evnum);

        System.out.print("REC::Particle\n");
        ev.read(parBank);
        parBank.show();
        //printPions(parBank);
        System.out.println("");

        /*
        System.out.print("MC::Particle\n");
        ev.read(genBank);
        //genBank.show();
        printPions(genBank);
        System.out.println("");

        System.out.print("MC::Lund\n");
        ev.read(lundBank);
        //lundBank.show();
        printPions(lundBank);
        System.out.println("");
        */

        System.out.print("REC::Calorimeter\n");
        ev.read(calBank);
        calBank.show();
        System.out.println("");

        System.out.print("REC::Track\n");
        ev.read(trkBank);
        trkBank.show();
        System.out.println("");

        System.out.print("REC::Traj\n");
        ev.read(trajBank);
        trajBank.show();
        System.out.println("");
      };
    };
  };

  public static void printPions(Bank b) {
    float px,py,pz;
    int pid;
    boolean found;


    for(int k=0; k<b.getRows(); k++) { 
      pid = b.getInt("pid",k);
      found = false;
      if(Math.abs(pid)==211) found=true; // pions
      if(pid==113 || Math.abs(pid)==213) found=true; // rhos
      if(Math.abs(pid)==130) found=true; // K-long
      if(Math.abs(pid)==310) found=true; // K-short
      if(Math.abs(pid)==311 || Math.abs(pid)==321) found=true; // K0,K+,K-


      if(found) {
        px = b.getFloat("px",k); 
        py = b.getFloat("py",k); 
        pz = b.getFloat("pz",k); 
        /*
        System.out.print(evnum+" "+pid+" "+
          Math.sqrt( Math.pow(px,2) + Math.pow(py,2) + 
                     Math.pow(pz,2) + Math.pow(0.139571,2) ) +" "+
          Math.sqrt( Math.pow(px,2) + Math.pow(py,2) ) +"\n"
        );
        */
        System.out.print(evnum+" "+pid+" "+px+" "+py+" "+pz+"\n");
      };
    };
  };
    
};
