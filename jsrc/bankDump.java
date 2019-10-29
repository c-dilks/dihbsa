import org.jlab.jnp.hipo4.io.*;
import org.jlab.jnp.hipo4.data.*;

import java.lang.Math;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class bankDump{

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
    Bank recBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
    Bank genBank = new Bank(reader.getSchemaFactory().getSchema("MC::Particle"));
    Bank lundBank = new Bank(reader.getSchemaFactory().getSchema("MC::Lund"));
    Bank runconfigBank = new Bank(reader.getSchemaFactory().getSchema("RUN::config"));

    reader.getEvent(ev,0); // init

    int evnum;
    int npart;

    while(reader.hasNext()==true) {

      reader.nextEvent(ev);

      ev.read(runconfigBank);
      evnum = runconfigBank.getInt("event",0);

      if(whichEv==0 || whichEv==evnum) {
        System.out.println("\n-----");
        System.out.println("EVNUM = "+evnum);

        System.out.print("REC::Particle\t");
        ev.read(recBank);
        recBank.show();
        printPions(recBank);
        System.out.println("");

        System.out.print("MC::Particle\t");
        ev.read(genBank);
        genBank.show();
        printPions(genBank);
        System.out.println("");

        System.out.print("MC::Lund\t");
        ev.read(lundBank);
        lundBank.show();
        printPions(lundBank);
        System.out.println("");
      };
    };
  };

  public static void printPions(Bank b) {
    Float px,py,pz;

    for(int k=0; k<b.getRows(); k++) { 
      if(b.getInt("pid",k)==211 || b.getInt("pid",k)==-211) {
        px = b.getFloat("px",k); 
        py = b.getFloat("py",k); 
        pz = b.getFloat("pz",k); 
        System.out.print(b.getInt("pid",k)+" "+
          Math.sqrt( Math.pow(px,2) + Math.pow(py,2) + 
                     Math.pow(pz,2) + Math.pow(0.139571,2) ) +" "+
          Math.sqrt( Math.pow(px,2) + Math.pow(py,2) ) +"\n"
        );
      };
    };
  };
    
};
