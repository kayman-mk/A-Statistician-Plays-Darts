package de.matthiaskay.dart;

public class Stats {
    // board dimensions in mm (measured from center of board)
    public static final double R1 = 6.35; // distance to double bullseye
    public static final double R2 = 15.9; // distance to single bullseye
    public static final double R3 = 99;   // distance to inside triple ring
    public static final double R4 = 107;  // distance to outside triple ring
    public static final double R5 = 162;  // distance to inside double ring
    public static final double R = 170;   // distance to outside double ring

    // ordering of the numbers
    public static final int[] d = new int[] {20,1,18,4,13,6,10,15,2,17,3,19,7,16,8,11,14,9,12,5};

    // index ordering of the numbers (unit-based)
    public static final int[] ii = new int[] {2,9,11,4,20,6,13,15,18,7,16,19,5,17,8,14,10,3,12,1};

    //
    /////////////////////////////////////////////////////////////
    // EM for the simple model
    ////////////////////////////////////////////////////////////
    //

    public static double simpleEM(int[] scores) {
        return simpleEM(scores,100,100);
    }

    // Parameters
    // scores: an array of scores
    // sInit: initial guess for the variance
    // numIter: the number of iterations to run EM algorithm
    //
    // Returns an estimate of the variance
    public static double simpleEM(int[] scores, double sInit, int numIter) {
        double s = sInit;

        for (int i=0; i<numIter; i++) {
            s = simpleStep(scores,s);
        }

        return s;
    }

    // Parameters
    // scores: an array of scores
    // s: current estimate of variance
    private static double simpleStep(int[] scores, double s) {
        double[] a = new double[6];
        double[] b = new double[6];

        // calculate the integral (pg.6 in paper) over each region
        // note: the regions are defined as follows
        // 0 - double bullseye
        // 1 - single bullseye
        // 2 - single
        // 3 - double ring
        // 4 - triple ring
        // 5 - off the board
        a[0] = integNumerator(s,0,R1);
        a[1] = integNumerator(s,R1,R2);
        a[2] = integNumerator(s,R2,R3)/20 + integNumerator(s,R4,R5)/20;
        a[3] = integNumerator(s,R5,R)/20;
        a[4] = integNumerator(s,R3,R4)/20;
        a[5] = integNumerator(s,R,-1);

        // calculate the integral over each region
        b[0] = integDenominator(s,0,R1);
        b[1] = integDenominator(s,R1,R2);
        b[2] = integDenominator(s,R2,R3)/20 + integDenominator(s,R4,R5)/20;
        b[3] = integDenominator(s,R5,R)/20;
        b[4] = integDenominator(s,R3,R4)/20;
        b[5] = integDenominator(s,R,-1);

        int n = scores.length;
        double e = 0;
        for (int i=0; i<n; i++) {
            e += computeExp(scores[i], a, b);
        }

        return e/(2*n);
    }

    // Parameters
    // x: the score we observe
    // a: the numerators we computed for each region
    // b: the denominators for each region
    private static double computeExp(int x, double[] a, double[] b) {
        // return the appropriate expectation, based on how the score can be achieved
        if (x==1 || x==5 || x==7 || x==11 || x==13 || x==17 || x==19) {
            // these numbers can only be achieved by hitting a single
            return a[2]/b[2];
        }
        else if (x==2 || x==4 || x==8 || x==10 || x==14 || x==16 || x==20) {
            // single or double
            return (a[2]+a[3])/(b[2]+b[3]);
        }
        else if (x==3 || x==9 || x==15) {
            // single or triple region
            return (a[2]+a[4])/(b[2]+b[4]);
        }

        else if (x==6 || x==12 || x==18) {
            // single, double, or triple
            return (a[2]+a[3]+a[4])/(b[2]+b[3]+b[4]);
        }
        else if (x==24 || x==30 || x==36) {
            // double or triple
            return (a[3]+a[4])/(b[3]+b[4]);
        }
        else if (x==22 || x==26 || x==28 || x==32 || x==34 || x==38 || x==40) {
            // double only
            return a[3]/b[3];
        }
        else if (x==21 || x==27 || x==33 || x==39 || x==42 || x==45 || x==48 ||
                x==51 || x==54 || x==57 || x==60) {
            // triple only
            return a[4]/b[4];
        }
        else if (x==25) {
            // single bullseye
            return a[1]/b[1];
        }
        else if (x==50) {
            // double bullseye
            return a[0]/b[0];
        }
        else {
            // outside the board
            return a[5]/b[5];
        }
    }

    // Parameters
    // s: variance
    // r1: lower bound of integration
    // r2: upper bound of integration
    private static double integNumerator(double s, double r1, double r2) {
        if (r2 == -1) // r2 is assumed to be infinity
            return (r1*r1+2*s)*Math.exp(-r1*r1/(2*s));
        else
            return (r1*r1+2*s)*Math.exp(-r1*r1/(2*s)) - (r2*r2+2*s)*Math.exp(-r2*r2/(2*s));
    }

    // Parameters
    // s: variance
    // r1: lower bound of integration
    // r2: upper bound of integration
    private static double integDenominator(double s, double r1, double r2) {
        if (r2 == -1) // r2 is assumed to be infinity
            return Math.exp(-r1*r1/(2*s));
        else
            return Math.exp(-r1*r1/(2*s)) - Math.exp(-r2*r2/(2*s));
    }

    //
    //////////////////////////////////////
    ///////////////////////
    // EM for the general model
    /////////////////////////////////////////////////////////////
    //

    // number of importance samples
    private static final int NUM_IMP_SAMPLES = 1000;

    // random number generator for importance samplce
    private static java.util.Random rand = new java.util.Random(0);

    // board areas used for importance sampling
    private static final double AR_1 = Math.PI*R1*R1;              // double bullseye
    private static final double AR_2 = Math.PI*(R2*R2 - R1*R1);    // single bullseye
    private static final double AR_3 = Math.PI*(R3*R3 - R2*R2)/20; // lower single
    private static final double AR_4 = Math.PI*(R5*R5 - R4*R4)/20; // upper single
    private static final double AR_5 = Math.PI*(R*R - R5*R5)/20;   // double
    private static final double AR_6 = Math.PI*(R4*R4 - R3*R3)/20; // triple

    public static double[] generalEM(int[] scores) {
        return generalEM(scores,100,100,0,150);
    }

    /**
     *
     * @param scores scores array
     * @param s1Init initial values for s1 (x variance)
     * @param s2Init initial values for s2 (y variance)
     * @param s3Init initial values for s3 (correlation)
     * @param numIter number of steps to run the EM algorithm
     *
     * @return an array of the estimates: {s1, s2, s3}
     */
    public static double[] generalEM(int[] scores, double s1Init, double s2Init,
                                     double s3Init, int numIter) {
        double s1 = s1Init;
        double s2 = s2Init;
        double s3 = s3Init;

        rand.setSeed(0);

        for (int i=0; i<numIter; i++) {
            double[] a = generalStep(scores, s1, s2, s3);
            s1 = a[0];
            s2 = a[1];
            s3 = a[2];
        }

        return new double[] {s1, s2, s3};
    }

    /**
     *
     * @param scores scores array
     * @param s1 current estimate of covariance matrix
     * @param s2 current estimate of covariance matrix
     * @param s3 current estimate of covariance matrix
     *
     * @return an array of updated estimates {s1, s2, s3}
     */
    public static double[] generalStep(int[] scores, double s1, double s2, double s3) {

        int n = scores.length;
        double[] a = new double[3];
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;

        for (int i=0; i<n; i++) {
            double[] b = simulateExp(scores[i], s1, s2, s3);
            a[0] += b[0];
            a[1] += b[1];
            a[2] += b[2];
        }

        return new double[] {a[0]/n, a[1]/n, a[2]/n};
    }

    // Parameters
    // x: the score we observe
    // s1, s2, s3: current estimate of covariance
    private static double[] simulateExp(int x, double s1, double s2, double s3) {
        double det = s1*s2 - s3*s3;
        double[] a = new double[3];
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;
        double big_w = 0;

        for (int i = 0; i< NUM_IMP_SAMPLES; i++) {
            double[] z = randomPt(x);
            double w = Math.exp(-(s2*z[0]*z[0] - 2*s3*z[0]*z[1] +
                    s1*z[1]*z[1])/(2*det));
            a[0] += w*z[0]*z[0];
            a[1] += w*z[1]*z[1];
            a[2] += w*z[0]*z[1];
            big_w += w;
        }

        a[0] /= big_w;
        a[1] /= big_w;
        a[2] /= big_w;
        return a;
    }

    // Parameters
    // x: the score we observe
    private static double[] randomPt(int x) {
        double u = rand.nextDouble();

        if (x==1 || x==5 || x==7 || x==11 || x==13 || x==17 || x==19) {
            // single only
            if (u <= AR_3 /(AR_3 + AR_4)) {
                return randomSlicePt(x,R2,R3);
            }
            else {
                return randomSlicePt(x,R4,R5);
            }
        }

        else if (x==2 || x==4 || x==8 || x==10 || x==14 || x==16 || x==20) {
            // single or double
            if (u <= AR_5 /(AR_5 + AR_3 + AR_4)) {
                return randomSlicePt(x/2,R5,R);
            }
            else if (u <= (AR_5 + AR_3)/(AR_5 + AR_3 + AR_4)) {
                return randomSlicePt(x,R2,R3);
            }
            else {
                return randomSlicePt(x,R4,R5);
            }
        }

        else if (x==3 || x==9 || x==15) {
            // single or triple
            if (u <= AR_6 /(AR_6 + AR_3 + AR_4)) {
                return randomSlicePt(x/3,R3,R4);
            }
            else if (u <= (AR_6 + AR_3)/(AR_6 + AR_3 + AR_4)) {
                return randomSlicePt(x,R2,R3);
            }
            else {
                return randomSlicePt(x,R4,R5);
            }
        }

        else if (x==6 || x==12 || x==18) {
            // single, double, or triple
            if (u <= AR_6 /(AR_6 + AR_5 + AR_3 + AR_4)) {
                return randomSlicePt(x/3,R3,R4);
            }
            else if (u <= (AR_6 + AR_5)/
                    (AR_6 + AR_5 + AR_3 + AR_4)) {
                return randomSlicePt(x/2,R5,R);
            }
            else if (u <= (AR_6 + AR_5 + AR_3)/
                    (AR_6 + AR_5 + AR_3 + AR_4)) {
                return randomSlicePt(x,R2,R3);
            }
            else {
                return randomSlicePt(x,R4,R5);
            }
        }

        else if (x==24 || x==30 || x==36) {
            // double or triple
            if (u <= AR_6 /(AR_6 + AR_5)) {
                return randomSlicePt(x/3,R3,R4);
            }
            else {
                return randomSlicePt(x/2,R5,R);
            }
        }

        else if (x==22 || x==26 || x==28 || x==32 || x==34 || x==38 || x==40) {
            // double only
            return randomSlicePt(x/2,R5,R);
        }

        else if (x==21 || x==27 || x==33 || x==39 || x==42 || x==45 || x==48 ||
                x==51 || x==54 || x==57 || x==60) {
            // triple only
            return randomSlicePt(x/3,R3,R4);
        }

        else if (x==25) {
            // single bullseye
            return randomCirclePt(R1,R2);
        }

        else if (x==50) {
            // double bullseye
            return randomCirclePt(0,R1);
        }

        else {
            // outside the board
            // a bit of a cheat, the second number should be infinity
            return randomCirclePt(R,10*R);
        }
    }

    // Parameters
    // x: the score we observe
    // r1, r2: lower and upper limits
    private static double[] randomSlicePt(int x, double r1, double r2) {
        int k = ii[x-1];
        double u = rand.nextDouble();
        double th = -2*Math.PI/40 + (k-1)*2*Math.PI/20 + 2*Math.PI/20*u;
        th = Math.PI/2 - th;
        double r = randomR(r1,r2);
        return new double[] {r*Math.cos(th), r*Math.sin(th)};
    }

    // Parameters:
    // r1, r2: lower and upper limits
    private static double[] randomCirclePt(double r1, double r2) {
        double u = rand.nextDouble();
        double th = 2*Math.PI*u;
        double r = randomR(r1,r2);
        return new double[] {r*Math.cos(th), r*Math.sin(th)};
    }

    // Parameters
    // r1, r2: lower and upper limits
    private static double randomR(double r1, double r2) {
        double u = rand.nextDouble();
        return Math.sqrt(r1*r1 + (r2*r2-r1*r1)*u);
    }

    //
    /////////////////////////////////////////////////////////////
    // Funtions to compute grid of expected scores
    /////////////////////////////////////////////////////////////
    //

    // powers of 2
    private static final int[] pow2 = new int[] {1,2,4,8,16,32,64,128,256,512,1024};

    // Parameters
    // s: variance
    // n: size of grid (assuming square, and a power of 2)
    //
    // Returns an n x n grid of expected scores
    public static double[][] computeExpScores(double s, int n) {
        // round up to the nearest power of 2
        int m = 1;
        for (int i=0; i<pow2.length; i++) {
            m = pow2[i];
            if (n <= pow2[i]) {
                break;
            }
        }

        // note: due to the design of all the FFT code (taken from Numerical
        // Recipes in C), all of our arrays must be unit-based!

        // number of millimeters per pixel
        double c = 2*R/n;

        // build the scores array
        double[][] big_s = new double[2*m+1][2*m+1];
        for (int i=1; i<=2*m; i++) {
            for (int j=1; j<=2*m; j++) {
                big_s[i][j] = score((i-m)*c,(j-m)*c);
            }
        }

        // helpful constants
        int n1 = (int)Math.floor(n/2.0);
        int n2 = (int)Math.ceil(n/2.0);

        // if the variance is 0, just return the scores array!
        if (s == 0) {
            double[][] E = new double[n][n];
            for (int i=m-n1+1; i<=m+n2; i++) {
                for (int j=m-n1+1; j<=m+n2; j++) {
                    E[i-m+n1-1][j-m+n1-1] = big_s[i][j];
                }
            }
            return E;
        }

        // build the Gaussian density arrays
        double[] g1 = new double[2*m+1];
        double[] g2 = new double[2*m+1];
        for (int i=1; i<=2*m; i++) {
            g1[i] = Math.exp(-(i-m)*(i-m)*c*c/(2*s))/Math.sqrt(2*Math.PI*s)*c;
            g2[i] = g1[i];
        }

        // compute the FT of g1 and g2
        double[] g1f = new double[4*m+1];
        double[] g2f = new double[4*m+1];
        twofft(g1,g2,g1f,g2f,2*m);

        // compute the FT of each row of S
        double[][] A = new double[2*m+1][4*m+1];
        for (int i=1; i<=2*m-1; i+=2) {
            twofftrow(big_s,i,i+1,A,2*m);
        }

        // multiply every row of A by g2f
        double re,im;
        for (int i=1; i<=2*m; i+=1) {
            for (int j=1; j<=4*m-1; j+=2) {
                re = A[i][j]*g2f[j] - A[i][j+1]*g2f[j+1];
                im = A[i][j]*g2f[j+1] + A[i][j+1]*g2f[j];
                A[i][j] = re;
                A[i][j+1] = im;
            }
        }

        // compute the inverse FT of every row of A
        for (int i=1; i<=2*m; i++) {
            four1row(A,i,2*m,-1);
        }

        // take the real part of each row, and switch around some columns
        double[][] AA = new double[2*m+1][2*m+1];
        for (int i=1; i<=2*m; i++) {
            for (int j=1; j<=m; j++) {
                AA[i][j+m] = A[i][2*j-1]/(2*m);
            }
            for (int j=m+1; j<=2*m; j++) {
                AA[i][j-m] = A[i][2*j-1]/(2*m);
            }
        }

        // compute the FT of every column of AA
        double[][] AAA = new double[4*m+1][2*m+1];
        for (int j=1; j<=2*m-1; j+=2) {
            twofftcol(AA,j,j+1,AAA,2*m);
        }

        // multiply every column of AAA by g1f
        for (int j=1; j<=2*m; j++) {
            for (int i=1; i<=4*m-1; i+=2) {
                re = AAA[i][j]*g1f[i] - AAA[i+1][j]*g1f[i+1];
                im = AAA[i][j]*g1f[i+1] + AAA[i+1][j]*g1f[i];
                AAA[i][j] = re;
                AAA[i+1][j] = im;
            }
        }

        // compute the inverse FT of every column of AAA
        for (int j=1; j<=2*m; j++) {
            four1col(AAA,j,2*m,-1);
        }

        // take the real part of each column, switch around some rows,
        // and strip away some part of the array on each side
        double[][] E = new double[n][n];
        for (int i=1; i<=n1; i++) {
            for (int j=m-n1+1; j<=m+n2; j++) {
                E[i+n2-1][j-m+n1-1] = AAA[2*i-1][j]/(2*m);
            }
        }
        for (int i=2*m-n2+1; i<=2*m; i++) {
            for (int j=m-n1+1; j<=m+n2; j++) {
                E[i-2*m+n2-1][j-m+n1-1] = AAA[2*i-1][j]/(2*m);
            }
        }

        // help the gc out
        big_s=null; g1=null; g2=null; g1f=null; g2f=null;
        A=null; AA=null; AAA=null;

        return E;
    }

    // Parameters
    // s1, s2, s3: x variance, y variance, covariance
    // n: size of grid (assuming square, and a power of 2)
    //
    // Returns an n x n grid of expected scores
    //
    // Note: as opposed to the simple version of this function:
    // computeExpScores(double s, int n)
    // this function is not really optimized for speed. The rationale
    // here is that it will take long enough to run the general EM
    // algorithm anyway.
    public static double[][] computeExpScores(double s1, double s2, double s3, int n) {
        // round up to the nearest power of 2
        int m = 1;
        for (int i=0; i<pow2.length; i++) {
            m = pow2[i];
            if (n <= pow2[i]) {
                break;
            }
        }

        // note: due to the design of all the FFT code (taken from Numerical
        // Recipes in C), all of our arrays must be unit-based!

        // another note: all of our matrices here are complex. the layout is that
        // there are twice as many columns as a corresponding real array representing
        // the same dimension, so that m[i][j] and m[i][j+1] represent corresponding
        // real and complex parts of the same number

        // number of millimeters per pixel
        double c = 2*R/n;

        // build the scores array
        double[][] S = new double[2*m+1][4*m+1];
        for (int i=1; i<=2*m; i++) {
            for (int j=1; j<=2*m; j++) {
                S[i][2*j-1] = score((i-m)*c,(j-m)*c);
                S[i][2*j] = 0;
            }
        }

        // helpful constants
        int n1 = (int)Math.floor(n/2.0);
        int n2 = (int)Math.ceil(n/2.0);

        // build the Gaussian density array
        double[][] G = new double[2*m+1][4*m+1];
        double det = s1*s2-s3*s3;
        for (int i=1; i<=2*m; i++) {
            for (int j=1; j<=2*m; j++) {
                G[i][2*j-1] = Math.exp(-(s2*(i-m)*(i-m) -
                        2*s3*(i-m)*(j-m) +
                        s1*(j-m)*(j-m))*c*c/(2*det))
                        /(2*Math.PI*Math.sqrt(det))*c*c;
                G[i][2*j] = 0;
            }
        }

        // compute the FT of each matrix (in-place)
        four2(S,2*m,1);
        four2(G,2*m,1);

        // multiply S and G together (store the result in S)
        double re,im;
        for (int i=1; i<=2*m; i+=1) {
            for (int j=1; j<=4*m-1; j+=2) {
                re = S[i][j]*G[i][j] - S[i][j+1]*G[i][j+1];
                im = S[i][j]*G[i][j+1] + S[i][j+1]*G[i][j];
                S[i][j] = re;
                S[i][j+1] = im;
            }
        }

        // compute the inverse FT of S (in-place)
        four2(S,2*m,-1);

        // switch around some rows and some columns, take the real part
        double[][] EE = new double[2*m+1][2*m+1];
        for (int i=1; i<=m; i++) {
            for (int j=1; j<=m; j++) {
                EE[i+m][j+m] = S[i][2*j-1]/(4*m*m);
            }
            for (int j=m+1; j<=2*m; j++) {
                EE[i+m][j-m] = S[i][2*j-1]/(4*m*m);
            }
        }
        for (int i=m+1; i<=2*m; i++) {
            for (int j=1; j<=m; j++) {
                EE[i-m][j+m] = S[i][2*j-1]/(4*m*m);
            }
            for (int j=m+1; j<=2*m; j++) {
                EE[i-m][j-m] = S[i][2*j-1]/(4*m*m);
            }
        }

        // trim away on each side
        double[][] E = new double[n][n];
        for (int i=m-n1; i<=m+n2-1; i++) {
            for (int j=m-n1; j<=m+n2-1; j++) {
                E[i-m+n1][j-m+n1] = EE[i][j];
            }
        }

        return E;
    }

    public static int score(double x, double y) {
        // compute the radius
        double r = Math.sqrt(x*x + y*y);

        // check if it's off the board (do this for speed)
        if (r > R) return 0;

        // check for a double bullseye
        if (r <= R1) return 50;

        // check for a single bullseye
        if (r <= R2) return 25;

        // get the angle
        double th = Math.atan2(y, x);
        double phi = (Math.PI/2+Math.PI/20-th) % (2*Math.PI);
        if (phi < 0) phi += 2*Math.PI;

        // now get the number
        int i = (int)(phi/(2*Math.PI)*20);
        if (i < 0) i = 0;
        if (i >= 19) i = 19;
        int n = d[i];

        // check for a single
        if (r <= R3) return n;

        // check for a triple
        if (r <= R4) return 3*n;

        // check for a single
        if (r <= R5) return n;

        //f we got here, it must be a double
        return 2*n;
    }

    public static String getRegion(double x, double y) {
        // compute the radius
        double r = Math.sqrt(x*x + y*y);

        // check if it's off the board (do this for speed)
        if (r > R) return "Off";

        // check for a double bullseye
        if (r <= R1) return "DB";

        // check for a single bullseye
        if (r <= R2) return "SB";

        // get the angle
        double th = Math.atan2(y, x);
        double phi = (Math.PI/2+Math.PI/20-th) % (2*Math.PI);
        if (phi < 0) phi += 2*Math.PI;

        // now get the number
        int i = (int)(phi/(2*Math.PI)*20);
        if (i < 0) i = 0;
        if (i >= 19) i = 19;
        int n = d[i];

        // check for a single
        if (r <= R3) return "S"+n;

        // check for a triple=
        if (r <= R4) return "T"+n;

        // check for a single
        if (r <= R5) return "S"+n;

        //f we got here, it must be a double
        return "D"+n;
    }

    public static String getRegion(int i, int j, int n1, int n2) {
        return getRegion((i+1-n1/2)*2*R/n1, (j+2-n2/2)*2*R/n2);
    }

    public static double[] getMaxAndArgmax(double[][] data) {
        double max = data[0][0];
        int ii=0, jj=0;

        for (int i=0; i<data.length; i++) {
            for (int j=0; j<data[0].length; j++) {
                if (data[i][j] > max) {
                    max = data[i][j];
                    ii = i;
                    jj = j;
                }
            }
        }

        return new double[] {max, ii, jj};
    }

    //
    /////////////////////////////////////////////////////////////
    // FFT Funtions
    // Adapted from Numerical Recipes in C
    /////////////////////////////////////////////////////////////
    //

    /**
     * Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1;
     * or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform, if
     * isign is input as -1.
     *
     * @param data is a complex array of length nn or, equivalently, a real array of length 2*nn
     * @param nn MUST be an integer power of 2 (this is not checked for!).
     */
    public static void four1(double[] data, int nn, int isign) {
        int n,mmax,m,j,istep,i;
        double wtemp,wr,wpr,wpi,wi,theta;
        double temp,tempr,tempi;

        n=nn << 1;
        j=1;
        for (i=1;i<n;i+=2) {
            if (j > i) {
                // SWAP(data[j],data[i]);
                temp=data[j];
                data[j]=data[i];
                data[i]=temp;
                // SWAP(data[j+1],data[i+1]);
                temp=data[j+1];
                data[j+1]=data[i+1];
                data[i+1]=temp;
            }
            m=n >>> 1;
            while (m >= 2 && j > m) {
                j -= m;
                m=m >>> 1;
            }
            j += m;
        }
        mmax=2;
        while (n > mmax) {
            istep=mmax << 1;
            theta=isign*(6.28318530717959/mmax);
            wtemp=Math.sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi=Math.sin(theta);
            wr=1.0;
            wi=0.0;
            for (m=1;m<mmax;m+=2) {
                for (i=m; i<=n;i+=istep) {
                    j=i+mmax;
                    tempr=wr*data[j]-wi*data[j+1];
                    tempi=wr*data[j+1]+wi*data[j];
                    data[j]=data[i]-tempr;
                    data[j+1]=data[i+1]-tempi;
                    data[i] += tempr;
                    data[i+1] += tempi;
                }
                wr=(wtemp=wr)*wpr-wi*wpi+wr;
                wi=wi*wpr+wtemp*wpi+wi;
            }
            mmax=istep;
        }
    }

    // four1 applied to the kth row of data (a matrix)
    private static void four1row(double[][] data, int k, int nn, int isign) {
        int n,mmax,m,j,istep,i;
        double wtemp,wr,wpr,wpi,wi,theta;
        double temp,tempr,tempi;

        n=nn << 1;
        j=1;
        for (i=1;i<n;i+=2) {
            if (j > i) {
                // SWAP(data[k][j],data[k][i]);
                temp=data[k][j];
                data[k][j]=data[k][i];
                data[k][i]=temp;
                // SWAP(data[k][j+1],data[k][i+1]);
                temp=data[k][j+1];
                data[k][j+1]=data[k][i+1];
                data[k][i+1]=temp;
            }
            m=n >>> 1;
            while (m >= 2 && j > m) {
                j -= m;
                m=m >>> 1;
            }
            j += m;
        }
        mmax=2;
        while (n > mmax) {
            istep=mmax << 1;
            theta=isign*(6.28318530717959/mmax);
            wtemp=Math.sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi=Math.sin(theta);
            wr=1.0;
            wi=0.0;
            for (m=1;m<mmax;m+=2) {
                for (i=m; i<=n;i+=istep) {
                    j=i+mmax;
                    tempr=wr*data[k][j]-wi*data[k][j+1];
                    tempi=wr*data[k][j+1]+wi*data[k][j];
                    data[k][j]=data[k][i]-tempr;
                    data[k][j+1]=data[k][i+1]-tempi;
                    data[k][i] += tempr;
                    data[k][i+1] += tempi;
                }
                wr=(wtemp=wr)*wpr-wi*wpi+wr;
                wi=wi*wpr+wtemp*wpi+wi;
            }
            mmax=istep;
        }
    }

    // four1 applied to the kth column of data (a matrix)
    private static void four1col(double[][] data, int k, int nn, int isign) {
        int n,mmax,m,j,istep,i;
        double wtemp,wr,wpr,wpi,wi,theta;
        double temp,tempr,tempi;

        n=nn << 1;
        j=1;
        for (i=1;i<n;i+=2) {
            if (j > i) {
                // SWAP(data[j][k],data[i][k]);
                temp=data[j][k];
                data[j][k]=data[i][k];
                data[i][k]=temp;
                // SWAP(data[j+1][k],data[i+1][k]);
                temp=data[j+1][k];
                data[j+1][k]=data[i+1][k];
                data[i+1][k]=temp;
            }
            m=n >>> 1;
            while (m >= 2 && j > m) {
                j -= m;
                m=m >>> 1;
            }
            j += m;
        }
        mmax=2;
        while (n > mmax) {
            istep=mmax << 1;
            theta=isign*(6.28318530717959/mmax);
            wtemp=Math.sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi=Math.sin(theta);
            wr=1.0;
            wi=0.0;
            for (m=1;m<mmax;m+=2) {
                for (i=m; i<=n;i+=istep) {
                    j=i+mmax;
                    tempr=wr*data[j][k]-wi*data[j+1][k];
                    tempi=wr*data[j+1][k]+wi*data[j][k];
                    data[j][k]=data[i][k]-tempr;
                    data[j+1][k]=data[i+1][k]-tempi;
                    data[i][k] += tempr;
                    data[i+1][k] += tempi;
                }
                wr=(wtemp=wr)*wpr-wi*wpi+wr;
                wi=wi*wpr+wtemp*wpi+wi;
            }
            mmax=istep;
        }
    }

    // Given two real input arrays data1[1..n] and data2[1..n], this routine calls four1
    // and returns two complex output arrays, fft1[1..2n] and fft2[1..2n], each of
    // complex length n (i.e., real length 2*n), which contain the discrete Fourier
    // transforms of the respective data arrays. n MUST be an integer power of 2.
    private static void twofft(double[] data1, double[] data2,
                               double[] fft1, double[] fft2, int n) {
        int nn3,nn2,jj,j;
        double rep,rem,aip,aim;

        nn3=1+(nn2=2+n+n);
        for (j=1,jj=2;j<=n;j++,jj+=2) {
            fft1[jj-1]=data1[j];
            fft1[jj]=data2[j];
        }
        four1(fft1,n,1);
        fft2[1]=fft1[2];
        fft1[2]=fft2[2]=0.0;
        for (j=3;j<=n+1;j+=2) {
            rep=0.5*(fft1[j]+fft1[nn2-j]);
            rem=0.5*(fft1[j]-fft1[nn2-j]);
            aip=0.5*(fft1[j+1]+fft1[nn3-j]);
            aim=0.5*(fft1[j+1]-fft1[nn3-j]);
            fft1[j]=rep;
            fft1[j+1]=aim;
            fft1[nn2-j]=rep;
            fft1[nn3-j] = -aim;
            fft2[j]=aip;
            fft2[j+1] = -rem;
            fft2[nn2-j]=aip;
            fft2[nn3-j]=rem;
        }
    }

    // twofft applied to the k1st and k2nd rows of data, results are stored
    // in the corresponding rows of fft
    private static void twofftrow(double[][] data, int k1, int k2,
                                  double[][] fft, int n) {
        int nn3,nn2,jj,j;
        double rep,rem,aip,aim;

        nn3=1+(nn2=2+n+n);
        for (j=1,jj=2;j<=n;j++,jj+=2) {
            fft[k1][jj-1]=data[k1][j];
            fft[k1][jj]=data[k2][j];
        }
        four1row(fft,k1,n,1);
        fft[k2][1]=fft[k1][2];
        fft[k1][2]=fft[k2][2]=0.0;
        for (j=3;j<=n+1;j+=2) {
            rep=0.5*(fft[k1][j]+fft[k1][nn2-j]);
            rem=0.5*(fft[k1][j]-fft[k1][nn2-j]);
            aip=0.5*(fft[k1][j+1]+fft[k1][nn3-j]);
            aim=0.5*(fft[k1][j+1]-fft[k1][nn3-j]);
            fft[k1][j]=rep;
            fft[k1][j+1]=aim;
            fft[k1][nn2-j]=rep;
            fft[k1][nn3-j] = -aim;
            fft[k2][j]=aip;
            fft[k2][j+1] = -rem;
            fft[k2][nn2-j]=aip;
            fft[k2][nn3-j]=rem;
        }
    }

    // twofft applied to the k1st and k2nd columns of data, results are stored
    // in the corresponding columns of fft
    private static void twofftcol(double[][] data, int k1, int k2,
                                  double[][] fft, int n) {
        int nn3,nn2,jj,j;
        double rep,rem,aip,aim;

        nn3=1+(nn2=2+n+n);
        for (j=1,jj=2;j<=n;j++,jj+=2) {
            fft[jj-1][k1]=data[j][k1];
            fft[jj][k1]=data[j][k2];
        }
        four1col(fft,k1,n,1);
        fft[1][k2]=fft[2][k1];
        fft[2][k1]=fft[2][k2]=0.0;
        for (j=3;j<=n+1;j+=2) {
            rep=0.5*(fft[j][k1]+fft[nn2-j][k1]);
            rem=0.5*(fft[j][k1]-fft[nn2-j][k1]);
            aip=0.5*(fft[j+1][k1]+fft[nn3-j][k1]);
            aim=0.5*(fft[j+1][k1]-fft[nn3-j][k1]);
            fft[j][k1]=rep;
            fft[j+1][k1]=aim;
            fft[nn2-j][k1]=rep;
            fft[nn3-j][k1] = -aim;
            fft[j][k2]=aip;
            fft[j+1][k2] = -rem;
            fft[nn2-j][k2]=aip;
            fft[nn3-j][k2]=rem;
        }
    }

    // Replaces data[1..n][1..2*nn] by its discrete 2D Fourier transform, if isign
    // is 1; or replaces data by nn times its inverse discrete Fourier transform, if
    // isign is -1. note: data is a complex matrix.
    // n MUST be an integer power of 2 (this is not checked for!).
    private static void four2(double[][] data, int nn, int isign) {
        // compute the 1D FT of every row
        for (int i=1; i<=nn; i++) {
            four1row(data,i,nn,isign);
        }

        // compute the 1D FT of every column
        double[] a = new double[2*nn+1];
        for (int j=1; j<=2*nn-1; j+=2) {
            // copy into the array a
            for (int i=1; i<=nn; i++) {
                a[2*i-1] = data[i][j];
                a[2*i] = data[i][j+1];
            }
            // compute the 1D FT of a
            four1(a,nn,isign);
            // copy back into the matrix data
            for (int i=1; i<=nn; i++) {
                data[i][j] = a[2*i-1];
                data[i][j+1] = a[2*i];
            }
        }
    }

    //
    /////////////////////////////////////////////////////////////
    // Functions for the DartsApplet class to call
    ////////////////////////////////////////////////////////////
    //

    public static double simpleEM(int[] scores, DartsApplet da) {
        return simpleEM(scores,100,100,da);
    }

    public static double simpleEM(int[] scores, double sInit, int numIter,
                                  DartsApplet da) {
        double s = sInit;

        int inc = numIter/10;

        for (int i=0; i<numIter; i++) {
            s = simpleStep(scores, s);
            if (i % inc == 0) {
                da.updateProgress(50*i/numIter);
            }
        }

        return s;
    }

    public static double[] generalEM(int[] scores, DartsApplet da) {
        return generalEM(scores,100,100,0,150,da);
    }

    public static double[] generalEM(int[] scores, double s1Init, double s2Init,
                                     double s3Init, int numIter, DartsApplet da) {
        double s1 = s1Init;
        double s2 = s2Init;
        double s3 = s3Init;

        rand.setSeed(0);

        int inc = numIter/10;

        for (int i=0; i<numIter; i++) {
            double[] A = generalStep(scores, s1, s2, s3);
            s1 = A[0];
            s2 = A[1];
            s3 = A[2];
            if (i % inc == 0) {
                da.updateProgress(50*i/numIter);
            }
        }

        return new double[] {s1, s2, s3};
    }

    public static double[][] computeExpScores(double s, int n, DartsApplet da) {
        // round up to the nearest power of 2
        int m = 1;
        for (int i=0; i<pow2.length; i++) {
            m = pow2[i];
            if (n <= pow2[i]) {
                break;
            }
        }

        // note: due to the design of all the FFT code (taken from Numerical
        // Recipes in C), all of our arrays must be unit-based!

        // number of millimeters per pixel
        double c = 2*R/n;

        // build the scores array
        double[][] S = new double[2*m+1][2*m+1];
        for (int i=1; i<=2*m; i++) {
            for (int j=1; j<=2*m; j++) {
                S[i][j] = score((i-m)*c,(j-m)*c);
            }
        }
        da.updateProgress(50+1/11*50);

        // helpful constants
        int n1 = (int)Math.floor(n/2.0);
        int n2 = (int)Math.ceil(n/2.0);

        // if the variance is 0, just return the scores array!
        if (s == 0) {
            double[][] E = new double[n][n];
            for (int i=m-n1+1; i<=m+n2; i++) {
                for (int j=m-n1+1; j<=m+n2; j++) {
                    E[i-m+n1-1][j-m+n1-1] = S[i][j];
                }
            }
            da.updateProgress(100);
            return E;
        }

        // build the Gaussian density arrays
        double[] g1 = new double[2*m+1];
        double[] g2 = new double[2*m+1];
        for (int i=1; i<=2*m; i++) {
            g1[i] = Math.exp(-(i-m)*(i-m)*c*c/(2*s))/Math.sqrt(2*Math.PI*s)*c;
            g2[i] = g1[i];
        }
        da.updateProgress(50+2/11*50);

        // compute the FT of g1 and g2
        double[] g1f = new double[4*m+1];
        double[] g2f = new double[4*m+1];
        twofft(g1,g2,g1f,g2f,2*m);
        da.updateProgress(50+3/11*50);

        // compute the FT of each row of S
        double[][] A = new double[2*m+1][4*m+1];
        for (int i=1; i<=2*m-1; i+=2) {
            twofftrow(S,i,i+1,A,2*m);
        }
        da.updateProgress(50+4/11*50);

        // multiply every row of A by g2f
        double re,im;
        for (int i=1; i<=2*m; i+=1) {
            for (int j=1; j<=4*m-1; j+=2) {
                re = A[i][j]*g2f[j] - A[i][j+1]*g2f[j+1];
                im = A[i][j]*g2f[j+1] + A[i][j+1]*g2f[j];
                A[i][j] = re;
                A[i][j+1] = im;
            }
        }
        da.updateProgress(50+5/11*50);

        // compute the inverse FT of every row of A
        for (int i=1; i<=2*m; i++) {
            four1row(A,i,2*m,-1);
        }
        da.updateProgress(50+6/11*50);

        // take the real part of each row, and switch around some columns
        double[][] AA = new double[2*m+1][2*m+1];
        for (int i=1; i<=2*m; i++) {
            for (int j=1; j<=m; j++) {
                AA[i][j+m] = A[i][2*j-1]/(2*m);
            }
            for (int j=m+1; j<=2*m; j++) {
                AA[i][j-m] = A[i][2*j-1]/(2*m);
            }
        }
        da.updateProgress(50+7/11*50);

        // compute the FT of every column of AA
        double[][] AAA = new double[4*m+1][2*m+1];
        for (int j=1; j<=2*m-1; j+=2) {
            twofftcol(AA,j,j+1,AAA,2*m);
        }
        da.updateProgress(50+8/11*50);

        // multiply every column of AAA by g1f
        for (int j=1; j<=2*m; j++) {
            for (int i=1; i<=4*m-1; i+=2) {
                re = AAA[i][j]*g1f[i] - AAA[i+1][j]*g1f[i+1];
                im = AAA[i][j]*g1f[i+1] + AAA[i+1][j]*g1f[i];
                AAA[i][j] = re;
                AAA[i+1][j] = im;
            }
        }
        da.updateProgress(50+9/11*50);

        // compute the inverse FT of every column of AAA
        for (int j=1; j<=2*m; j++) {
            four1col(AAA,j,2*m,-1);
        }
        da.updateProgress(50+10/11*50);

        // take the real part of each column, switch around some rows,
        // and strip away some part of the array on each side
        double[][] E = new double[n][n];
        for (int i=1; i<=n1; i++) {
            for (int j=m-n1+1; j<=m+n2; j++) {
                E[i+n2-1][j-m+n1-1] = AAA[2*i-1][j]/(2*m);
            }
        }
        for (int i=2*m-n2+1; i<=2*m; i++) {
            for (int j=m-n1+1; j<=m+n2; j++) {
                E[i-2*m+n2-1][j-m+n1-1] = AAA[2*i-1][j]/(2*m);
            }
        }
        da.updateProgress(100);

        // help the gc out
        S=null; g1=null; g2=null; g1f=null; g2f=null;
        A=null; AA=null; AAA=null;

        return E;
    }

    public static double[][] computeExpScores(double s1, double s2, double s3,
                                              int n, DartsApplet da) {
        // round up to the nearest power of 2
        int m = 1;
        for (int i=0; i<pow2.length; i++) {
            m = pow2[i];
            if (n <= pow2[i]) {
                break;
            }
        }

        // note: due to the design of all the FFT code (taken from Numerical
        // Recipes in C), all of our arrays must be unit-based!

        // another note: all of our matrices here are complex. the layout is that
        // there are twice as many columns as a corresponding real array representing
        // the same dimension, so that m[i][j] and m[i][j+1] represent corresponding
        // real and complex parts of the same number

        // number of millimeters per pixel
        double c = 2*R/n;

        // build the scores array
        double[][] S = new double[2*m+1][4*m+1];
        for (int i=1; i<=2*m; i++) {
            for (int j=1; j<=2*m; j++) {
                S[i][2*j-1] = score((i-m)*c,(j-m)*c);
                S[i][2*j] = 0;
            }
        }
        da.updateProgress(50+1/7*50);

        // helpful constants
        int n1 = (int)Math.floor(n/2.0);
        int n2 = (int)Math.ceil(n/2.0);

        // build the Gaussian density array
        double[][] G = new double[2*m+1][4*m+1];
        double det = s1*s2-s3*s3;
        for (int i=1; i<=2*m; i++) {
            for (int j=1; j<=2*m; j++) {
                G[i][2*j-1] = Math.exp(-(s2*(i-m)*(i-m) -
                        2*s3*(i-m)*(j-m) +
                        s1*(j-m)*(j-m))*c*c/(2*det))
                        /(2*Math.PI*Math.sqrt(det))*c*c;
                G[i][2*j] = 0;
            }
        }
        da.updateProgress(50+2/7*50);

        // compute the FT of each matrix (in-place)
        four2(S,2*m,1);
        da.updateProgress(50+3/7*50);

        four2(G,2*m,1);
        da.updateProgress(50+4/7*50);

        // multiply S and G together (store the result in S)
        double re,im;
        for (int i=1; i<=2*m; i+=1) {
            for (int j=1; j<=4*m-1; j+=2) {
                re = S[i][j]*G[i][j] - S[i][j+1]*G[i][j+1];
                im = S[i][j]*G[i][j+1] + S[i][j+1]*G[i][j];
                S[i][j] = re;
                S[i][j+1] = im;
            }
        }
        da.updateProgress(50+5/7*50);

        // compute the inverse FT of S (in-place)
        four2(S,2*m,-1);
        da.updateProgress(50+6/7*50);

        // switch around some rows and some columns, take the real part
        double[][] EE = new double[2*m+1][2*m+1];
        for (int i=1; i<=m; i++) {
            for (int j=1; j<=m; j++) {
                EE[i+m][j+m] = S[i][2*j-1]/(4*m*m);
            }
            for (int j=m+1; j<=2*m; j++) {
                EE[i+m][j-m] = S[i][2*j-1]/(4*m*m);
            }
        }
        for (int i=m+1; i<=2*m; i++) {
            for (int j=1; j<=m; j++) {
                EE[i-m][j+m] = S[i][2*j-1]/(4*m*m);
            }
            for (int j=m+1; j<=2*m; j++) {
                EE[i-m][j-m] = S[i][2*j-1]/(4*m*m);
            }
        }

        // trim away on each side
        double[][] E = new double[n][n];
        for (int i=m-n1; i<=m+n2-1; i++) {
            for (int j=m-n1; j<=m+n2-1; j++) {
                E[i-m+n1][j-m+n1] = EE[i][j];
            }
        }
        da.updateProgress(100);

        // help the gc out
        S=null; G=null; EE=null;

        return E;
    }

    public static void main(String[] args) {

        int[] andyScores = new int[] {12,16,19,3,17,1,25,19,17,50,18,1,3,17,2,2,13,18,16,2,25,5,5,1,5,4,17,25,25,50,3,7,17,17,3,3,3,7,11,10,25,1,19,15,4,1,5,12,17,16,50,20,20,20,25,50,2,17,3,20,20,20,5,1,18,15,2,3,25,12,9,3,3,19,16,20,5,5,1,4,15,16,5,20,16,2,25,6,12,25,11,25,7,2,5,19,17,17,2,12,12,9,14,25,25,20,1,20,50,2};
        int[] ryanScores = new int[] {20,20,5,19,11,25,15,2,7,1,12,1,5,19,3,17,17,2,10,16,17,19,5,20,8,0,17,17,2,15,2,25,60,19,17,2,4,1,5,15,16,7,2,13,5,0,7,19,14,2,11,18,20,19,17,15,11,5,5,4,9,20,1,17,2,0,19,3,2,13,4,4,7,7,3,17,14,12,19,17,17,17,13,15,14,19,3,6,13,18,9,20,9,12,19,18,0,12,17,17,3,5,8,19,3,17,10,15,8,51};
        int[] jonScores = new int[] {19,1,18,12,12,14,12,5,16,19,19,15,1,6,14,19,3,9,4,17,16,7,16,28,20,20,2,19,3,0,16,5,15,17,17,19,9,18,4,2,19,19,5,1,18,15,17,3,2,2,2,10,3,20,2,15,19,18,20,5,12,12,18,17,11,14,25,5,5,1,19,19,20,25,15,10,16,8,15,3,19,9,4,4,9,17,2,10,11,5,15,9,2,17,17,3,8,18,19,15,4,13,13,4,14,17,2,16,9,5,10,15,17,3,9,17,19,16,16,4,33,8,19,17,4,5,50,9,2,17,3,7,20,20,20,1,19,2,5,18,8,7,3,2,15,15,40,17,3,19};
        System.out.println(Math.sqrt(Stats.simpleEM(andyScores)));
        double[] p = Stats.generalEM(andyScores);
        System.out.println(Math.sqrt(p[0]));
        System.out.println(Math.sqrt(p[1]));
        System.out.println(p[2]/Math.sqrt(p[0]*p[1]));
    }
}
