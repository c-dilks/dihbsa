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

    Event ev = new Event();
    Bank recBank = new Bank(reader.getSchemaFactory().getSchema("REC::Particle"));
    Bank genBank = new Bank(reader.getSchemaFactory().getSchema("MC::Particle"));
    Bank lundBank = new Bank(reader.getSchemaFactory().getSchema("MC::Lund"));
    Bank runconfigBank = new Bank(reader.getSchemaFactory().getSchema("RUN::config"));

    reader.getEvent(ev,0); // init

    int pid = -211;

    while(reader.hasNext()==true) {

      reader.nextEvent(ev);

      ev.read(runconfigBank);
      System.out.println("\n-----");
      System.out.println("EVNUM = "+runconfigBank.getInt("event",0));

      System.out.print("REC::Particle\t");
      ev.read(recBank);
      recBank.show();
      //for(int k=0; k<recBank.getRows(); k++) { if(recBank.getInt("pid",k)==pid) System.out.print(recBank.getFloat("pz",k)+" "); };
      System.out.println("");

      System.out.print("MC::Particle\t");
      ev.read(genBank);
      genBank.show();
      //for(int k=0; k<genBank.getRows(); k++) { if(genBank.getInt("pid",k)==pid) System.out.print(genBank.getFloat("pz",k)+" "); };
      System.out.println("");

      System.out.print("MC::Lund\t");
      ev.read(lundBank);
      lundBank.show();
      //for(int k=0; k<lundBank.getRows(); k++) { if(lundBank.getInt("pid",k)==pid) System.out.print(lundBank.getFloat("pz",k)+" "); };
      System.out.println("");
    };
  };
};
