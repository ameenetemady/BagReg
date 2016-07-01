package proteinidentification;

import java.util.ArrayList;

public class RAFeature {

    static final int FEATURE_NUMBER = 5;
    private int matchedPeptide;
    private int uniquePeptide;
    private int matchedSpectrum;
    private int uniqueSpectrum;
    private double maxScore;
    private double averageScore;
    private double criterion;
    private String criterionName;

    public RAFeature() {
        criterion = 0;
        criterionName = "matchedPeptide";
    }

    public RAFeature(String proteinName) {
        matchedPeptide = RAFileManager.getMatchedPeptideNumber(proteinName);
        uniquePeptide = RAFileManager.getUniquePeptideNumber(proteinName);
        matchedSpectrum = RAFileManager.getMatchedSpectrumNumber(proteinName);
        uniqueSpectrum = RAFileManager.getUniqueSpectrumNumber(proteinName);
        maxScore = RAFileManager.getMaxScore(proteinName);
        averageScore = RAFileManager.getAverageScore(proteinName);

        criterion = 0;
    }

    public double[] getFeatureValuesBesidesCriterion() {
        double[] featureSet = new double[FEATURE_NUMBER-1];
        ArrayList<String> names = getFeatureNames();
        int j = 0;
        for (int i = 0; i < names.size(); i++) {
            if (!names.get(i).equals(criterionName)) {
                featureSet[j++] = getFeatureValueByName(names.get(i));
            }
        }
        return featureSet;
    }

    public static ArrayList<String> getFeatureNames() {
        ArrayList<String> names = new ArrayList<String>();
        names.add("matchedPeptide");
        names.add("uniquePeptide");
        names.add("matchedSpectrum");
//        names.add("uniqueSpectrum");
        names.add("maxScore");
        names.add("averageScore");
        return names;
    }

    public void setCriterion(String featureName) {
        criterionName = featureName;
        if (featureName.equals("matchedPeptide")) {
            criterion = matchedPeptide;
        } else if (featureName.equals("uniquePeptide")) {
            criterion = uniquePeptide;
        } else if (featureName.equals("matchedSpectrum")) {
            criterion = matchedSpectrum;
        } else if (featureName.equals("uniqueSpectrum")) {
            criterion = uniqueSpectrum;
        } else if (featureName.equals("maxScore")) {
            criterion = maxScore;
        } else {
            criterion = averageScore;
        }
    }

    public double getCriterion() {
        return criterion;
    }

    private double getFeatureValueByName(String featureName) {
        if (featureName.equals("matchedPeptide")) {
            return 1.0 * matchedPeptide;
        } else if (featureName.equals("uniquePeptide")) {
            return 1.0 * uniquePeptide;
        } else if (featureName.equals("matchedSpectrum")) {
            return 1.0 * matchedSpectrum;
        } else if (featureName.equals("uniqueSpectrum")) {
            return 1.0 * uniqueSpectrum;
        } else if (featureName.equals("maxScore")) {
            return maxScore;
        } else {
            return averageScore;
        }
    }

    public void setMatchedPeptide(int mp) {
        matchedPeptide = mp;
    }

    public void setUniquePeptide(int up) {
        uniquePeptide = up;
    }

    public void setMatchedSpectrum(int ms) {
        matchedSpectrum = ms;
    }

    public void setUniqueSpectrum(int us) {
        uniqueSpectrum = us;
    }

    public void setMaxScore(double ms) {
        maxScore = ms;
    }

    public void setAverageScore(double as) {
        averageScore = as;
    }

    public void printFeatureValues() {
            System.out.println(this.matchedPeptide + "(mp) " + this.uniquePeptide + "(up) " + this.matchedSpectrum + "(ms) "
                    + this.uniqueSpectrum + "(us) " + this.maxScore + "(max) " + this.averageScore + "(ave)");
    }
}
