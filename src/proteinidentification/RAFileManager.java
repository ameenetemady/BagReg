package proteinidentification;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

public class RAFileManager {

    public static String sampleName;
    public static String inputPath;

    private static String getProteinName(String s) {
        String proteinName;
        if (s.contains("sp|") || s.contains("gi|")) {
            proteinName = s.split("\\|")[1].trim();
        } else if (s.contains("[Contaminant]")) {
            proteinName = s.substring(13).trim();
        } else if (s.contains("|")) {
            proteinName = s.split("\\|")[0].trim();
        } else {
            proteinName = s.trim();
        }
        return proteinName;
    }

    public static HashSet<String> getProteinSet() {
        HashSet<String> set = new HashSet<String>();
        try {
            File file = new File(inputPath);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String string;
            while ((string = br.readLine()) != null) {
                String[] cut = string.split(" +");
                String proteinName = getProteinName(cut[2].trim());
                set.add(proteinName);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace(System.err);
            System.exit(0);
        } finally {
            return set;
        }
    }

    public static int getMatchedPeptideNumber(String name) {
        HashSet<String> peptideSet = new HashSet<String>();
        try {
            File file = new File(inputPath);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String string;
            while ((string = br.readLine()) != null) {
                String[] cut = string.split(" +");
                String proteinName, peptide;
                peptide = cut[1].trim();
                proteinName = getProteinName(cut[2].trim());
                if (proteinName.equals(name)) {
                    peptideSet.add(peptide);
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace(System.err);
            System.exit(0);
        } finally {
            return peptideSet.size();
        }
    }

    public static int getUniquePeptideNumber(String name) {
        HashSet<String> peptideSet = new HashSet<String>();
        Hashtable<String, String> matchedProtein = new Hashtable<String, String>();
        try {
            File file = new File(inputPath);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String string;
            while ((string = br.readLine()) != null) {
                String[] cut = string.split(" +");
                String proteinName, peptide;
                peptide = cut[1].trim();
                proteinName = getProteinName(cut[2].trim());
                if (proteinName.equals(name)) {
                    peptideSet.add(peptide);
                }
                String proteinList;
                if (matchedProtein.containsKey(peptide)) {
                    proteinList = matchedProtein.get(peptide);
                    if (!proteinList.contains(proteinName)) {
                        proteinList += " " + proteinName;
                    }
                } else {
                    proteinList = proteinName;
                }
                matchedProtein.put(peptide, proteinList);
            }
            br.close();
            Iterator<String> iterator = peptideSet.iterator();
            while (iterator.hasNext()) {
                String peptide = iterator.next();
                String proteinList = matchedProtein.get(peptide);
                if (!proteinList.equals(name)) {
                    iterator.remove();
                }
            }
        } catch (Exception e) {
            e.printStackTrace(System.err);
            System.exit(0);
        } finally {
            return peptideSet.size();
        }
    }

    public static int getMatchedSpectrumNumber(String name) {
        HashSet<String> spectrumSet = new HashSet<String>();
        try {
            File file = new File(inputPath);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String string;
            while ((string = br.readLine()) != null) {
                String[] cut = string.split(" +");
                String proteinName, spectrum;
                spectrum = cut[0].trim();
                proteinName = getProteinName(cut[2].trim());
                if (proteinName.equals(name)) {
                    if (spectrumSet.contains(spectrum)) {
//                        System.out.println(spectrum + "\t" + proteinName);
                    }
                    spectrumSet.add(spectrum);
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace(System.err);
            System.exit(0);
        } finally {
            return spectrumSet.size();
        }
    }

    public static int getUniqueSpectrumNumber(String name) {
        HashSet<String> totalSpectrumSet = new HashSet<String>();
        HashSet<String> duplicateSpectrumSet = new HashSet<String>();
        HashSet<String> peptideSet = new HashSet<String>();
//        Hashtable<String, String> matchedProtein = new Hashtable<String, String>();
        try {
            File file = new File(inputPath);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String string;
            while ((string = br.readLine()) != null) {
                String[] cut = string.split(" +");
                String proteinName, spectrum, peptide;
                spectrum = cut[0].trim();
                peptide = cut[1].trim();
                proteinName = getProteinName(cut[2].trim());
                if (proteinName.equals(name)) {
                    totalSpectrumSet.add(spectrum);
                    if (peptideSet.contains(peptide)) {
                        duplicateSpectrumSet.add(spectrum);
                    } else {
                        peptideSet.add(peptide);
                    }
                }
//                String proteinList;
//                if (matchedProtein.containsKey(spectrum)) {
//                    proteinList = matchedProtein.get(spectrum);
//                    proteinList += " ";
//                } else {
//                    proteinList = new String();
//                }
//                proteinList += proteinName;
//                matchedProtein.put(spectrum, proteinList);
            }
            br.close();
//            Iterator<String> iterator = spectrumSet.iterator();
//            while (iterator.hasNext()) {
//                String spectrum = iterator.next();
//                String proteinList = matchedProtein.get(spectrum);
//                if (!proteinList.equals(name)) {
//                    iterator.remove();
//                }
//            }
        } catch (Exception e) {
            e.printStackTrace(System.err);
            System.exit(0);
        } finally {
            return totalSpectrumSet.size() - duplicateSpectrumSet.size();
        }
    }

    public static double getMaxScore(String name) {
        double maxScore = 0;
        try {
            File file = new File(inputPath);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String string;
            while ((string = br.readLine()) != null) {
                String[] cut = string.split(" +");
                String proteinName;
                double score;
                proteinName = getProteinName(cut[2].trim());
                score = Double.parseDouble(cut[3].trim());
                if (proteinName.equals(name)) {
                    if (score > maxScore) {
                        maxScore = score;
                    }
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace(System.err);
            System.exit(0);
        } finally {
            return maxScore;
        }
    }

    public static double getAverageScore(String name) {
        Hashtable<String, Double> peptideScore = new Hashtable<String, Double>();
        double averageScore = 0;
        try {
            File file = new File(inputPath);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String string;
            while ((string = br.readLine()) != null) {
                String[] cut = string.split(" +");
                String proteinName, peptide;
                Double score;
                peptide = cut[1].trim();
                proteinName = getProteinName(cut[2].trim());
                score = Double.valueOf(cut[3].trim());
                if (proteinName.equals(name)) {
                    if (peptideScore.containsKey(peptide)) {
                        Double temp = peptideScore.get(peptide);
                        if (temp.compareTo(score) < 0) {
                            peptideScore.put(peptide, score);
                        }
                    } else {
                        peptideScore.put(peptide, score);
                    }
                }
            }
            br.close();
            double sumScore = 0;
            Enumeration<String> e = peptideScore.keys();
            while (e.hasMoreElements()) {
                String peptide = e.nextElement();
                sumScore += peptideScore.get(peptide).doubleValue();
            }
            averageScore = sumScore / peptideScore.size();
        } catch (Exception e) {
            e.printStackTrace(System.err);
            System.exit(0);
        } finally {
            return averageScore;
        }
    }

    public static RAFeature getFeature(String name) {
        RAFeature feature = new RAFeature();
        HashSet<String> peptideSet = new HashSet<String>();
        Hashtable<String, HashSet<String>> matchedProtein = new Hashtable<String, HashSet<String>>();
        HashSet<String> totalSpectrumSet = new HashSet<String>();
        Hashtable<String, Double> peptideScore = new Hashtable<String, Double>();
        Hashtable<String, Integer> peptideStatistic = new Hashtable<String, Integer>();
        double averageScore = 0, maxScore = 0;
        int lineNum = 0;
        try {
            File file = new File(inputPath);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String string;
            while ((string = br.readLine()) != null) {
                lineNum++;
                String[] cut = string.split(" +");
                String proteinName, peptide, spectrum;
                Double score;
                spectrum = cut[0].trim();
                peptide = cut[1].trim();
                proteinName = getProteinName(cut[2].trim());
                score = Double.valueOf(cut[3].trim());
                if (proteinName.equals(name)) {
                    if (!totalSpectrumSet.contains(spectrum)) {
                        Integer integer = 0;
                        if (peptideStatistic.containsKey(peptide)) {
                            integer = peptideStatistic.get(peptide);
                        }
                        integer++;
                        peptideStatistic.put(peptide, integer);
                    }
                    totalSpectrumSet.add(spectrum);
                    peptideSet.add(peptide);
                    if (score.doubleValue() > maxScore) {
                        maxScore = score;
                    }
                    if (peptideScore.containsKey(peptide)) {
                        Double temp = peptideScore.get(peptide);
                        if (temp.compareTo(score) < 0) {
                            peptideScore.put(peptide, score);
                        }
                    } else {
                        peptideScore.put(peptide, score);
                    }
                }
                HashSet<String> proteinSet=null;
                if (matchedProtein.containsKey(peptide)) {
                    proteinSet = matchedProtein.get(peptide);
                } else {
                    proteinSet = new HashSet<String>();
                }
                proteinSet.add(proteinName);
                matchedProtein.put(peptide, proteinSet);
            }
//            System.out.println("Line number:"+lineNum);
            br.close();
            double sumScore = 0;
            Enumeration<String> e = peptideScore.keys();
            while (e.hasMoreElements()) {
                String peptide = e.nextElement();
                sumScore += peptideScore.get(peptide).doubleValue();
            }
            averageScore = sumScore / peptideScore.size();
            
            feature.setMatchedPeptide(peptideSet.size());
            feature.setMatchedSpectrum(totalSpectrumSet.size());
            feature.setMaxScore(maxScore);
            feature.setAverageScore(averageScore);

            int count = 0;
            e=peptideStatistic.keys();
            while(e.hasMoreElements()){
                String peptide = e.nextElement();
                Integer integer = peptideStatistic.get(peptide);
                if(integer.equals(new Integer(1))){
                    count++;
                }
            }
            feature.setUniqueSpectrum(count);
            count = 0;
            e = matchedProtein.keys();
            while(e.hasMoreElements()){
                String peptide = e.nextElement();
                HashSet<String> proteinSet = matchedProtein.get(peptide);
                if(proteinSet.size() == 1 && proteinSet.contains(name)){
                    count++;
                }
            }
            feature.setUniquePeptide(count);
        } catch (Exception e) {
            e.printStackTrace(System.err);
            System.exit(0);
        } finally {
            return feature;
        }
    }

    public static void writeResultToFile(String output, ArrayList<Map.Entry<String, Double>> list) {
        PrintWriter writer = null;
        try {
            File writeFile = new File(output);
            if (writeFile.exists() == false) {
                writeFile.createNewFile();
                writeFile = new File(output);
            }
//            System.out.println("********************RESULT************************");
            writer = new PrintWriter(new FileOutputStream(writeFile));
            for (int i = 0; i < list.size(); i++) {
                writer.println(list.get(i).getKey() + "," + list.get(i).getValue());
//                System.out.println(list.get(i).getKey()+"\t"+list.get(i).getValue());
                writer.flush();
            }
        } catch (Exception x) {
            x.printStackTrace(System.err);
        } finally {
            if (writer != null) {
                writer.close();
            }
        }
    }

    public static void writeROCToFile(String output, List<Roc> ROC) {
        PrintWriter pw = null;
        File writefile;
        try {
            writefile = new File(output);
            if (writefile.exists() == false) {
                writefile.createNewFile();
                writefile = new File(output);
            }
            pw = new PrintWriter(new FileOutputStream(writefile));
            int i = 0;
            for(i = 0; i < ROC.size(); i++){
                pw.println(ROC.get(i).fpr+","+ROC.get(i).tpr);
                pw.flush();
            }
        } catch (Exception x) {
            x.printStackTrace(System.err);
        } finally {
            if (null != pw) {
                pw.close();
            }
        }
    }

    public static Hashtable<String, Double> getExistingResult(String input) {
        Hashtable<String, Double> result = new Hashtable<String, Double>();
        try {
            File file = new File(input);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String string;
            while ((string = br.readLine()) != null) {
                String[] cut = string.split(",");
                result.put(cut[0].trim(), Double.valueOf(cut[1].trim()));
//                System.out.println(cut[0].trim()+"\t"+result.get(cut[0].trim()));
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace(System.err);
            System.exit(0);
        } finally {
            return result;
        }
    }

    public static HashSet<String> getReferenceSet(String input) {
        HashSet<String> referenceSet = new HashSet<String>();
        try {
            File file = new File(input);
            BufferedReader br = new BufferedReader(new FileReader(file));
            String string;
            while ((string = br.readLine()) != null) {
                referenceSet.add(string.trim());
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace(System.err);
            System.exit(0);
        } finally {
            return referenceSet;
        }
    }
}
