package de.unikoeln.chemie.nmr.io;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import org.jcamp.parser.JCAMPException;
import org.jcamp.parser.JCAMPReader;
import org.jcamp.spectrum.ArrayData;
import org.jcamp.spectrum.IAssignmentTarget;
import org.jcamp.spectrum.IDataArray1D;
import org.jcamp.spectrum.IOrderedDataArray1D;
import org.jcamp.spectrum.NMRSpectrum;
import org.jcamp.spectrum.OrderedArrayData;
import org.jcamp.spectrum.Peak;
import org.jcamp.spectrum.Peak1D;
import org.jcamp.spectrum.assignments.AtomReference;
import org.jcamp.units.CommonUnit;
import org.jcamp.units.Unit;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;

import de.unikoeln.chemie.nmr.data.NmreData;

public class NmredataReader {
	BufferedReader input = null;
	
	Map<String,String> spectra=new HashMap<String,String>();
	Map<String,Peak> signals=new HashMap<String,Peak>();
	Map<String,IAssignmentTarget[]> assignments=new HashMap<String,IAssignmentTarget[]>();
	
	public NmredataReader(Reader in){
        input = new BufferedReader(in);
	}
	
	public NmredataReader(InputStream in){
		this(new InputStreamReader(in));
	}
	
	public NmreData read() throws Exception{
		NmreData data=new NmreData();
		IteratingSDFReader mdlreader=new IteratingSDFReader(input, DefaultChemObjectBuilder.getInstance());
		IAtomContainer ac = mdlreader.next();
		data.setMolecule(ac);
		mdlreader.close();
		String signalblock=null;
		for(Object key : ac.getProperties().keySet()){
			if(((String)key).startsWith("NMREDATA_1D")){
				spectra.put(((String)key).substring(9),(String)ac.getProperties().get(key));
			}else if(((String)key).equals("NMREDATA_ASSIGNMENT")){
				signalblock=(String)ac.getProperties().get(key);
			}
		}
		if(signalblock!=null)
			analyzeSignals(data, signalblock);
		else
			throw new Exception("There is no NMREDATA_ASSIGNMENT block in this file - required");
		analyzeSpectra(data);
		return data;
	}
	

	private void analyzeSignals(NmreData data, String signalblock) throws Exception {
		StringTokenizer st=new StringTokenizer(signalblock,"\n\r");
		while(st.hasMoreTokens()){
			String line = st.nextToken();
			StringTokenizer st2 = new StringTokenizer(line,",");
			String label=st2.nextToken();
			double shift = Double.parseDouble(st2.nextToken().trim());
			Peak peak=new Peak1D(shift,0);
			List<AtomReference> atoms = new ArrayList<AtomReference>();
			while(st2.hasMoreTokens()){
				String atom = st2.nextToken();
				if(atom.indexOf("H")>-1){
					int atomid=Integer.parseInt(atom.trim().substring(1))-1;
					if(atomid>=data.getMolecule().getAtomCount())
						throw new Exception("Atom "+atomid+" specified in MREDATA_ASSIGNMENT block, but only "+data.getMolecule().getAtomCount()+" atoms are in Molecule");
                    for(int k=0;k<data.getMolecule().getConnectedAtomsCount(data.getMolecule().getAtom(atomid));k++){
                        if(data.getMolecule().getConnectedAtomsList(data.getMolecule().getAtom(atomid)).get(k).getSymbol().equals("H")){
                        	atomid=data.getMolecule().getAtomNumber(data.getMolecule().getConnectedAtomsList(data.getMolecule().getAtom(atomid)).get(k));
                        	break;
                        }
                    }
					atoms.add(new AtomReference(null, atomid));
				}else{
					int atomid=Integer.parseInt(atom.trim())-1;
					atoms.add(new AtomReference(null, atomid));
				}
				IAssignmentTarget[] assigns = new IAssignmentTarget[atoms.size()];
				for (int i=0;i<atoms.size();i++)
					assigns[i] = (IAssignmentTarget) atoms.get(i);
			}
			signals.put(label, peak);
		}
	}

	private void analyzeSpectra(NmreData data) throws JCAMPException {
		for(String spectrum : spectra.keySet()){
			analyze1DSpectrum(spectra.get(spectrum), spectrum.substring(spectrum.indexOf("_")+1), data);
		}
		
	}

	private void analyze1DSpectrum(String spectrumblock, String nucleus, NmreData data) throws JCAMPException {
		StringTokenizer st=new StringTokenizer(spectrumblock,"\n\r");
		List<Peak> peaks=new ArrayList<>();
		double freq=Double.NaN;
		while(st.hasMoreTokens()){
			String line = st.nextToken();
			if(line.startsWith("Larmor=")){
				freq=Double.parseDouble(line.substring(7));
			}else{
				StringTokenizer st2 = new StringTokenizer(line,",");
				Peak peak=null;
				String multiplicity;
				while(st2.hasMoreTokens()){
					String token=st2.nextToken().trim();
					if(token.indexOf("=")==-1){
						//TODO use information
					}else if(token.startsWith("L")){
						peak=signals.get(token.substring(2).trim());
					}else if(token.startsWith("S")){
						multiplicity=token.substring(2).trim();
					}
				}
				//TODO multiplicity
				peaks.add(new Peak1D(peak.getPosition()[0],0));
			}
		}
        NMRSpectrum spectrum = null;
        Unit xUnit =  CommonUnit.hertz;
        Unit yUnit = CommonUnit.intensity;
        double reference = 0;
        Peak1D[] peaks1d = new Peak1D[peaks.size()];
        int i=0;
        for(Peak peak : peaks){
        	peaks1d[i]=(Peak1D)peak;
        	i++;
        }
        double[][] xy = peakTableToPeakSpectrum(peaks1d);
        IOrderedDataArray1D x = new OrderedArrayData(xy[0], xUnit);
        IDataArray1D y = new ArrayData(xy[1], yUnit);
        spectrum = new NMRSpectrum(x, y, nucleus, freq, reference, false, JCAMPReader.RELAXED);
        spectrum.setPeakTable(peaks1d);
        //spectrum.setAssignments((Assignment[]) tables[2]);
        data.addSpectrum(spectrum);
	}

	
	/**
	 * create peak spectrum from peak table.
	 * adds all intensities belonging to the same x-position up
	 * @param peaks Peak1D[]
	 * @return double[][] array of {x,  y}
	 */
	protected static double[][] peakTableToPeakSpectrum(Peak1D[] peaks)
	    throws JCAMPException {
	    int n = peaks.length;
	    if (n == 0)
	        throw new JCAMPException("empty peak table");
	    Arrays.sort(peaks);    
	    ArrayList<Double> px = new ArrayList<>(n);
	    ArrayList<Double> py = new ArrayList<>(n);
	    double x0 = peaks[0].getPosition()[0];
	    double y0 = peaks[0].getHeight();
	    for (int i = 1; i < n; i++) {
	        double x = peaks[i].getPosition()[0];
	        double y = peaks[i].getHeight();
	        if (x - x0 > Double.MIN_VALUE) {
	            px.add(new Double(x0));
	            py.add(new Double(y0));
	            x0 = x;
	            y0 = y;
	        } else {
	            y0 += y;
	        }
	    }
	    px.add(new Double(x0));
	    py.add(new Double(y0));
	    double[][] xy = new double[2][px.size()];
	    for (int i = 0; i < px.size(); i++) {
	        xy[0][i] = ((Double) px.get(i)).doubleValue();
	        xy[1][i] = ((Double) py.get(i)).doubleValue();
	    }
	    return xy;
	}
}
