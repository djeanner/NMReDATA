package de.unikoeln.chemie.nmr.test;

import java.io.IOException;
import java.io.InputStream;

import org.openscience.cdk.exception.CDKException;

import de.unikoeln.chemie.nmr.data.NmreData;
import de.unikoeln.chemie.nmr.io.NmredataReader;
import junit.framework.Assert;
import junit.framework.TestCase;

public class NmredataReaderTest  extends TestCase{
	
	public void testSimple() throws Exception, IOException{
		String filename = "testdata/example.nmredata.sdf";
        InputStream ins = this.getClass().getClassLoader().getResourceAsStream(filename);
        NmredataReader reader = new NmredataReader(ins);
        NmreData data = reader.read();
        Assert.assertEquals(32, data.getMolecule().getAtomCount());
        Assert.assertEquals(36, data.getMolecule().getBondCount());
        Assert.assertEquals(1, data.getSpectra().size());
        Assert.assertEquals(7, data.getSpectra().get(0).getPeakTable().length);
	}
	
}
