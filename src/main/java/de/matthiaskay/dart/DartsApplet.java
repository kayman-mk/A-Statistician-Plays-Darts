package de.matthiaskay.dart;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.util.*;
import javax.swing.colorchooser.*;

public class DartsApplet implements MouseMotionListener, ChangeListener {
    private JTextArea textArea;
    private JButton createHeatMapButton;
    private JLabel legendLabel;
    private JButton changeColorsButton;
    private DartsHeatMap heatMap;
    private JFrame changeColorsFrame;
    private JTextField regionField, expScoreField;
    private JCheckBox showBoardCheckBox, showArgmaxCheckBox;
    private JPanel gridPanel;
    private JTextField maxExpScoreField, aimForField;
    private JLabel sigmaXLabel, sigmaYLabel, rhoLabel;
    private JTextField sigmaXField, sigmaYField, rhoField;
    private JLabel numScoresLabel;
    private JRadioButton simpleButton, generalButton;
    private JProgressBar progressBar;
    private JLabel progressLabel;

    private boolean showBoard = true, showOptimum = true;
    private boolean simpleModel = true;

    public void init(JFrame window) {
        /////////////////////////////////
        // the west panel
        /////////////////////////////////

        // the instructions label
        JPanel instructionsPanel = new JPanel();
        instructionsPanel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Instructions"),
                BorderFactory.createEmptyBorder(5,5,0,5)));
        JLabel instructionsLabel = new JLabel("<html>" +
                "1. Throw 50 or so darts, aim-<br>&nbsp;&nbsp;&nbsp; " +
                "ing at the double bullseye.<br><br>" +
                "2. Enter the scores of each<br>&nbsp;&nbsp;&nbsp;throw below, " +
                "separated<br>&nbsp;&nbsp;&nbsp;by commas.<br><br>" +
                "3. Click \"Create heat map!\"<br>&nbsp;&nbsp;&nbsp;to see a " +
                "personalized<br>&nbsp;&nbsp;&nbsp;heat map." +
                "</html>");
        instructionsPanel.add(instructionsLabel);

        // text area to enter the scores
        textArea = new JTextArea();
        textArea.setLineWrap(true);
        JScrollPane scrollPane = new JScrollPane(textArea);
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
        scrollPane.setPreferredSize(new Dimension(200,300));
        scrollPane.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Scores"),
                BorderFactory.createEmptyBorder(5,5,5,5)));

        // create heat map button
        createHeatMapButton = new JButton("Create heat map!");
        createHeatMapButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                int r = createHeatMap();
                if (r > 0) {
                    String message = failMessages[r];
                    if (r == INVALID_NUMBERS_FAIL) {
                        message += invalidNumbers;
                    }
                    JOptionPane.showMessageDialog(window,
                            message, "Error", JOptionPane.ERROR_MESSAGE);
                }
            }
        });
        createHeatMapButton.setAlignmentX(JComponent.CENTER_ALIGNMENT);

        // put it together
        JPanel westPanel = new JPanel();
        westPanel.setPreferredSize(new Dimension(200,500));
        westPanel.setLayout(new BoxLayout(westPanel, BoxLayout.Y_AXIS));
        westPanel.add(instructionsPanel);
        westPanel.add(scrollPane);
        westPanel.add(Box.createVerticalStrut(5));
        westPanel.add(createHeatMapButton);
        westPanel.add(Box.createVerticalStrut(10));
        window.add(westPanel, BorderLayout.WEST);


        /////////////////////////////////
        // the main panel
        /////////////////////////////////

        // legend
        JPanel legendPanel = new JPanel();
        legendPanel.setLayout(new BoxLayout(legendPanel, BoxLayout.Y_AXIS));
        legendLabel = new JLabel();
        legendLabel.setBorder(BorderFactory.createLineBorder(Color.BLACK));
        legendLabel.setAlignmentX(JComponent.CENTER_ALIGNMENT);
        legendPanel.add(legendLabel);
        JLabel lowHighLabel = new JLabel("Low                 High");
        lowHighLabel.setAlignmentX(JComponent.CENTER_ALIGNMENT);
        legendPanel.add(lowHighLabel);

        // change colors
        changeColorsButton = new JButton("Change");
        changeColorsButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                changeColorsFrame.setVisible(true);
            }
        });
        buildChangeColorsFrame();
        JPanel changeColorsPanel = new JPanel();
        changeColorsPanel.setLayout(new BoxLayout(changeColorsPanel, BoxLayout.Y_AXIS));
        changeColorsPanel.add(changeColorsButton);
        changeColorsPanel.add(new JLabel(" "));

        // heat map
        heatMap = new DartsHeatMap();
        heatMap.setPreferredSize(new Dimension((int)(2*Stats.R),(int)(2*Stats.R)));
        heatMap.setBorder(BorderFactory.createLineBorder(Color.BLACK));
        // pay attention to the mouse
        heatMap.addMouseMotionListener(this);
        // a bit of a hack to make sure that the heat map size doesn't change
        JPanel heatMapPanel = new JPanel();
        heatMapPanel.add(heatMap);

        // set the color scheme
        setColorScheme(CLASSIC);

        JPanel colorSchemePanel = new JPanel();
        colorSchemePanel.add(new JLabel("<html>Color legend:<br>&nbsp;</html>"));
        colorSchemePanel.add(legendPanel);
        colorSchemePanel.add(Box.createHorizontalStrut(1));
        colorSchemePanel.add(changeColorsPanel);

        // mouse info panel
        JPanel mouseInfoPanel = new JPanel();
        mouseInfoPanel.add(new JLabel("Region:"));
        regionField = new JTextField(3);
        regionField.setEditable(false);
        mouseInfoPanel.add(regionField);
        mouseInfoPanel.add(Box.createHorizontalStrut(10));
        mouseInfoPanel.add(new JLabel("Expected score:"));
        expScoreField = new JTextField(4);
        expScoreField.setEditable(false);
        mouseInfoPanel.add(expScoreField);

        // options panel
        JPanel optionsPanel = new JPanel();
        showBoardCheckBox = new JCheckBox();
        showBoardCheckBox.setSelected(showBoard);
        showBoardCheckBox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                showBoard = !showBoard;
                heatMap.repaint();
            }
        });
        optionsPanel.add(showBoardCheckBox);
        optionsPanel.add(new JLabel("Show dartboard"));
        optionsPanel.add(Box.createHorizontalStrut(10));
        showArgmaxCheckBox = new JCheckBox();
        showArgmaxCheckBox.setSelected(showOptimum);
        showArgmaxCheckBox.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                showOptimum = !showOptimum;
                if (E != null) {
                    // smart repaint (for speed!)
                    // assume 1 mm per pixel
                    heatMap.repaint(0,imax-3,(int)(2*Stats.R)-1-jmax-3,6,6);
                }
            }
        });
        optionsPanel.add(showArgmaxCheckBox);
        optionsPanel.add(new JLabel("Show optimum place to aim"));

        // put it all together
        JPanel mainPanel = new JPanel();
        mainPanel.setPreferredSize(new Dimension(400,500));
        mainPanel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Heat map"),
                BorderFactory.createEmptyBorder(5,5,5,5)));
        mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));
        mainPanel.add(colorSchemePanel);
        mainPanel.add(heatMapPanel);
        mainPanel.add(mouseInfoPanel);
        mainPanel.add(Box.createVerticalGlue());
        mainPanel.add(optionsPanel);
        window.add(mainPanel, BorderLayout.CENTER);

        /////////////////////////////////
        // the east panel
        /////////////////////////////////

        // the help panel
        JPanel helpPanel = new JPanel();
        helpPanel.setPreferredSize(new Dimension(200,180));
        helpPanel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Help"),
                BorderFactory.createEmptyBorder(5,5,0,5)));
        helpPanel.add(new JLabel("<html>" +
                "First look at the color legend<br>to see which colors " +
                "corres-<br>pond to high expected scores.<br>Then look at the heat " +
                "map<br> to see which aiming locations<br>are assigned these " +
                "colors.<br>" +
                "(Note: the \u03c3 parameters are<br>measured " +
                "in mm; for refer-<br>ence, the dartboard's " +
                "radius<br>is 170 mm.)" +
                "</html>"));

        // the stats panel
        JPanel statsPanel = new JPanel();
        statsPanel.setPreferredSize(new Dimension(200,240));
        statsPanel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Stats"),
                BorderFactory.createEmptyBorder(5,5,0,5)));

        // the stats
        maxExpScoreField = new JTextField(4);
        maxExpScoreField.setEditable(false);
        maxExpScoreField.setForeground(Color.RED);
        aimForField = new JTextField(4);
        aimForField.setEditable(false);
        aimForField.setForeground(Color.RED);
        sigmaXLabel = new JLabel("\u03c3 = ");
        sigmaXField = new JTextField(4);
        sigmaXField.setEditable(false);
        sigmaXField.setForeground(Color.RED);
        sigmaYLabel = new JLabel("\u03c3_y =  ");
        sigmaYField = new JTextField(4);
        sigmaYField.setEditable(false);
        sigmaYField.setForeground(Color.RED);
        rhoLabel = new JLabel("\u03c1 = ");
        rhoField = new JTextField(4);
        rhoField.setEditable(false);
        rhoField.setForeground(Color.RED);

        gridPanel = new JPanel();
        gridPanel.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();
        c.anchor = GridBagConstraints.LINE_START;
        c.gridx = 0; c.gridy = 0;
        gridPanel.add(new JLabel("Max expected score: "),c);
        c.gridx = 1; c.gridy = 0;
        gridPanel.add(maxExpScoreField,c);
        c.gridx = 0; c.gridy = 1;
        gridPanel.add(new JLabel("Aim for: "),c);
        c.gridx = 1; c.gridy = 1;
        gridPanel.add(aimForField,c);
        c.gridx = 0; c.gridy = 2;
        gridPanel.add(sigmaXLabel,c);
        c.gridx = 1; c.gridy = 2;
        gridPanel.add(sigmaXField,c);

        // number of scores label
        numScoresLabel = new JLabel("<html>&nbsp;<br>&nbsp;</html>");
        c.insets = new Insets(7,0,0,0);
        c.gridx = 0; c.gridy = 7; c.gridwidth=2;
        gridPanel.add(numScoresLabel,c);

        // simple or general button
        simpleButton = new JRadioButton("Simple model");
        simpleButton.setSelected(true);
        simpleButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (!simpleModel) {
                    simpleModel = true;
                    sigmaXLabel.setText("\u03c3 = ");
                    // clear all the fields
                    maxExpScoreField.setText("");
                    aimForField.setText("");
                    sigmaXField.setText("");
                    sigmaYField.setText("");
                    rhoField.setText("");
                    // rearrange the grid panel
                    gridPanel.remove(sigmaYLabel);
                    gridPanel.remove(sigmaYField);
                    gridPanel.remove(rhoLabel);
                    gridPanel.remove(rhoField);
                    // clear the num scores label
                    numScoresLabel.setText("");
                }
            }
        });
        generalButton = new JRadioButton("General model");
        generalButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (simpleModel) {
                    simpleModel = false;
                    sigmaXLabel.setText("\u03c3_x = ");
                    // clear all the fields
                    maxExpScoreField.setText("");
                    aimForField.setText("");
                    sigmaXField.setText("");
                    sigmaYField.setText("");
                    rhoField.setText("");
                    // rearrange the grid panel
                    GridBagConstraints c = new GridBagConstraints();
                    c.anchor = GridBagConstraints.LINE_START;
                    c.gridx = 0; c.gridy = 3;
                    gridPanel.add(sigmaYLabel,c);
                    c.gridx = 1; c.gridy = 3;
                    gridPanel.add(sigmaYField,c);
                    c.gridx = 0; c.gridy = 4;
                    gridPanel.add(rhoLabel,c);
                    c.gridx = 1; c.gridy = 4;
                    gridPanel.add(rhoField,c);
                    // clear the num scores label
                    numScoresLabel.setText("");
                }
            }
        });
        ButtonGroup buttonGroup = new ButtonGroup();
        buttonGroup.add(simpleButton);
        buttonGroup.add(generalButton);

        c.insets = new Insets(7,0,0,0);
        c.gridx = 0; c.gridy = 8; c.gridwidth=2;
        gridPanel.add(simpleButton,c);

        c.insets = new Insets(0,0,0,0);
        c.gridx = 0; c.gridy = 9; c.gridwidth=2;
        gridPanel.add(generalButton,c);

        statsPanel.add(gridPanel);

        // the progress panel
        JPanel progressPanel = new JPanel();
        progressPanel.setPreferredSize(new Dimension(200,70));
        progressPanel.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Progress"),
                BorderFactory.createEmptyBorder(5,5,5,5)));
        progressPanel.setLayout(new BoxLayout(progressPanel, BoxLayout.Y_AXIS));

        // the progress label and bar
        progressLabel = new JLabel("<html>&nbsp;<br>&nbsp;</html>");
        progressPanel.add(progressLabel);

        progressBar = new JProgressBar(0,112);
        progressBar.setValue(0);
        progressBar.setStringPainted(true);
        JPanel progressBarPanel = new JPanel();
        progressBarPanel.add(progressBar);
        progressPanel.add(progressBarPanel);

        // put it together
        JPanel eastPanel = new JPanel();
        eastPanel.setPreferredSize(new Dimension(200,500));
        eastPanel.setLayout(new BoxLayout(eastPanel, BoxLayout.Y_AXIS));
        eastPanel.add(helpPanel);
        eastPanel.add(statsPanel);
        eastPanel.add(progressPanel);
        window.add(eastPanel, BorderLayout.EAST);
    }

    // pay attention to the mouse over the heat map
    // (we have to implement MouseMotionListener and add this as a separate function,
    // instead of the usual way, because otherwise this won't run on a JVM <= 1.5)
    public void mouseMoved(MouseEvent e) {
        if (E != null) {
            // assume 1 mm per pixel
            int i = e.getX();
            int j = (int)(2*Stats.R)-1-e.getY();
            int n = E.length;
            regionField.setText(Stats.getRegion(i,j,n,n));
            if (0 <= i && i < n && 0 <= j && j < n) {
                expScoreField.setText(Double.toString(roundTo2(E[i][j])));
            }
            else {
                expScoreField.setText("-----");
            }
        }
    }
    public void mouseDragged(MouseEvent e) {}

    private static double roundTo2(double d) {
        return Math.round(100*d)/100.0;
    }

    // the frame to change the color scheme
    JRadioButton classicButton, grayscaleButton, rainbowButton, customButton;
    JLabel[] customLabels;
    int numCustomColors;
    JColorChooser chooser;

    private void buildChangeColorsFrame() {
        changeColorsFrame = new JFrame("Change colors");

        JPanel choosePanel = new JPanel();
        choosePanel.add(new JLabel("Choose a color scheme:"));

        // classic
        classicButton = new JRadioButton("Classic");
        classicButton.setSelected(true);
        JPanel classicPreview = new JPanel();
        JLabel[] classicLabels = new JLabel[CLASSIC.length];
        for (int i=0; i<CLASSIC.length; i++) {
            classicLabels[i] = new JLabel("    ");
            classicLabels[i].setBorder(BorderFactory.createLineBorder(Color.BLACK));
            classicLabels[i].setBackground(CLASSIC[i]);
            classicLabels[i].setOpaque(true);
            classicPreview.add(classicLabels[i]);
        }

        // grayscale
        grayscaleButton = new JRadioButton("Grayscale");
        JPanel grayscalePreview = new JPanel();
        JLabel[] grayscaleLabels = new JLabel[GRAYSCALE.length];
        for (int i=0; i<GRAYSCALE.length; i++) {
            grayscaleLabels[i] = new JLabel("    ");
            grayscaleLabels[i].setBorder(BorderFactory.createLineBorder(Color.BLACK));
            grayscaleLabels[i].setBackground(GRAYSCALE[i]);
            grayscaleLabels[i].setOpaque(true);
            grayscalePreview.add(grayscaleLabels[i]);
        }

        // rainbow
        rainbowButton = new JRadioButton("Rainbow");
        JPanel rainbowPreview = new JPanel();
        JLabel[] rainbowLabels = new JLabel[RAINBOW.length];
        for (int i=0; i<RAINBOW.length; i++) {
            rainbowLabels[i] = new JLabel("    ");
            rainbowLabels[i].setBorder(BorderFactory.createLineBorder(Color.BLACK));
            rainbowLabels[i].setBackground(RAINBOW[i]);
            rainbowLabels[i].setOpaque(true);
            rainbowPreview.add(rainbowLabels[i]);
        }

        // custom
        customButton = new JRadioButton("Custom");
        JPanel customPreview = new JPanel();
        customLabels = new JLabel[10];
        for (int i=0; i<10; i++) {
            customLabels[i] = new JLabel("    ");
            customLabels[i].setBorder(BorderFactory.createLineBorder(Color.BLACK));
            customPreview.add(customLabels[i]);
        }
        customPreview.add(Box.createHorizontalStrut(3));
        numCustomColors = 0;
        JButton clearButton = new JButton("Clear");
        clearButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                for (int i=0; i<numCustomColors; i++) {
                    customLabels[i].setBackground(null);
                    customLabels[i].setOpaque(false);
                }
                numCustomColors = 0;
                customButton.setSelected(true);
            }
        });

        // button group
        ButtonGroup buttonGroup = new ButtonGroup();
        buttonGroup.add(classicButton);
        buttonGroup.add(grayscaleButton);
        buttonGroup.add(rainbowButton);
        buttonGroup.add(customButton);

        // put it all together
        JPanel radioPanel = new JPanel();
        radioPanel.setLayout(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();
        c.anchor = GridBagConstraints.LINE_START;
        c.gridx = 0; c.gridy = 0;
        radioPanel.add(classicButton,c);
        c.gridx = 1; c.gridy = 0;
        radioPanel.add(classicPreview,c);
        c.gridx = 0; c.gridy = 1;
        radioPanel.add(grayscaleButton,c);
        c.gridx = 1; c.gridy = 1;
        radioPanel.add(grayscalePreview,c);
        c.gridx = 0; c.gridy = 2;
        radioPanel.add(rainbowButton,c);
        c.gridx = 1; c.gridy = 2;
        radioPanel.add(rainbowPreview,c);
        c.gridx = 0; c.gridy = 3;
        radioPanel.add(customButton,c);
        c.gridx = 1; c.gridy = 3;
        radioPanel.add(customPreview,c);
        c.gridx = 2; c.gridy = 3;
        radioPanel.add(clearButton,c);

        // chooser
        chooser = new JColorChooser();
        chooser.setPreviewPanel(new JPanel());
        AbstractColorChooserPanel[] panels = chooser.getChooserPanels();
        for (int i=0; i<panels.length; i++) {
            String name = panels[i].getClass().getName();
            if (name.equals("javax.swing.colorchooser.DefaultRGBChooserPanel")) {
                chooser.removeChooserPanel(panels[i]);
            }
            else if (name.equals("javax.swing.colorchooser.DefaultHSBChooserPanel")) {
                chooser.removeChooserPanel(panels[i]);
            }
        }
        chooser.getSelectionModel().addChangeListener(this);

        // OK and cancel
        JPanel okCancelPanel = new JPanel();
        JButton okButton = new JButton("OK");
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (classicButton.isSelected()) {
                    setColorScheme(CLASSIC);
                }
                else if (grayscaleButton.isSelected()) {
                    setColorScheme(GRAYSCALE);
                }
                else if (rainbowButton.isSelected()) {
                    setColorScheme(RAINBOW);
                }
                else {
                    // if there were no custom colors, don't do anything
                    if (numCustomColors == 0) return;

                    Color[] customScheme = new Color[numCustomColors];
                    for (int i=0; i<numCustomColors; i++) {
                        customScheme[i] = customLabels[i].getBackground();
                    }
                    setColorScheme(customScheme);
                }
                heatMap.clear();
                heatMap.repaint();
                changeColorsFrame.setVisible(false);
            }
        });
        JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                closeChangeColorsFrame();
            }
        });
        okCancelPanel.add(okButton);
        okCancelPanel.add(cancelButton);

        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));
        mainPanel.setPreferredSize(new Dimension(450,320));
        mainPanel.add(choosePanel);
        mainPanel.add(radioPanel);
        mainPanel.add(chooser);
        mainPanel.add(okCancelPanel);
        changeColorsFrame.add(mainPanel, BorderLayout.CENTER);

        changeColorsFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        changeColorsFrame.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                closeChangeColorsFrame();
            }
        });
        changeColorsFrame.pack();
        changeColorsFrame.setResizable(false);
        changeColorsFrame.setVisible(false);
    }

    // for the color chooser
    public void stateChanged(ChangeEvent e) {
        if (numCustomColors < 10) {
            customLabels[numCustomColors].setBackground(chooser.getColor());
            customLabels[numCustomColors].setOpaque(true);
            numCustomColors++;
        }
        customButton.setSelected(true);
    }

    private void closeChangeColorsFrame() {
        if (heatMap.colorScheme == CLASSIC) {
            classicButton.setSelected(true);
        }
        else if (heatMap.colorScheme == GRAYSCALE) {
            grayscaleButton.setSelected(true);
        }
        else if (heatMap.colorScheme == RAINBOW) {
            rainbowButton.setSelected(true);
        }
        else {
            customButton.setSelected(true);
            for (int i=0; i<10; i++) {
                if (i < heatMap.colorScheme.length) {
                    customLabels[i].setBackground(heatMap.colorScheme[i]);
                    customLabels[i].setOpaque(true);
                }
                else {
                    customLabels[i].setBackground(null);
                    customLabels[i].setOpaque(false);
                }
            }
            numCustomColors = Math.min(heatMap.colorScheme.length,10);
        }
        changeColorsFrame.setVisible(false);
    }

    // heat map colors schemes
    private static final Color[] CLASSIC = new Color[]
            {new Color(255,0,0), new Color(255,113,0), new Color(255,255,0), new Color(255,255,255)};
    private static final Color[] GRAYSCALE = new Color[]
            {Color.BLACK, Color.WHITE};
    private static final Color[] RAINBOW = new Color[]
            {new Color(255,0,0), new Color(255,219,0), new Color(73,255,0), new Color(0,255,146),
                    new Color(0,146,255), new Color(73,0,255), new Color(255,0,219)};

    // the custom heat map component
    private class DartsHeatMap extends JPanel {
        public DartsHeatMap() {
            super();
            setBackground(Color.WHITE);
            setDoubleBuffered(true);
        }

        public Color[] colorScheme;
        private Color[] heatColors;
        private int numColors;

        public void setColorScheme(Color[] colors) {
            int k = colors.length;
            float[] r = new float[k];
            float[] g = new float[k];
            float[] b = new float[k];
            for (int j=0; j<k; j++) {
                r[j] = colors[j].getRed()/255.0f;
                g[j] = colors[j].getGreen()/255.0f;
                b[j] = colors[j].getBlue()/255.0f;
            }
            colorScheme = colors;
            numColors = Math.max(k*20,50);
            float m = numColors-1;
            float mm = k-1;

            heatColors = new Color[numColors];
            for (int i=0; i<numColors; i++) {
                float x = i/m*mm;
                int j = (int)(x);
                if (j <= k-2) {
                    heatColors[i] = new Color(r[j]*(1-x+j) + r[j+1]*(x-j),
                            g[j]*(1-x+j) + g[j+1]*(x-j),
                            b[j]*(1-x+j) + b[j+1]*(x-j));
                }
                else {
                    heatColors[i] = new Color(r[j],g[j],b[j]);
                }
            }
        }

        private Color getHeatColor(double t, double min, double max) {
            int i = (int)Math.round((t-min)/(max-min)*(numColors-1));
            return heatColors[i];
        }

        public BufferedImage makeLegend(int w, int h) {
            BufferedImage legend = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
            Graphics2D g2 = legend.createGraphics();

            for (int i=0; i<w; i++) {
                g2.setColor(getHeatColor(i,0,w-1));
                g2.fillRect(i,0,1,h);
            }

            return legend;
        }

        private BufferedImage bufferedImage;

        public void paintComponent(Graphics g) {
            super.paintComponent(g);

            int w = getWidth();
            int h = getHeight();
            Graphics2D g2 = (Graphics2D)g;

            // turn on antialiasing
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                    RenderingHints.VALUE_ANTIALIAS_ON);

            // clear the panel
            g2.setColor(Color.WHITE);
            g2.fillRect(0,0,w,h);

            // we'd like to make this image as few times as possible
            if (bufferedImage == null) {
                createBufferedImage(DartsApplet.this.E);
            }

            if (bufferedImage != null) {
                // draw the heat map
                g2.drawImage(bufferedImage,0,0,w,h,0,0,
                        bufferedImage.getWidth(),bufferedImage.getHeight(),null);

                // draw the board if we should
                if (DartsApplet.this.showBoard) {
                    drawBoard(g2,w,h);
                }

                // draw the optimum aiming spot if we should
                if (DartsApplet.this.showOptimum) {
                    drawOptimum(g2,w,h);
                }
            }
        }

        private void createBufferedImage(double[][] E) {
            if (E == null || E.length == 0) {
                bufferedImage = null;
                return;
            }

            int n = E.length;
            bufferedImage = new BufferedImage(n, n, BufferedImage.TYPE_INT_RGB);
            Graphics2D g2 = bufferedImage.createGraphics();

            // find the min and max values of E
            double min = E[0][0];
            double max = E[0][0];
            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    if (E[i][j] < min) min = E[i][j];
                    if (E[i][j] > max) max = E[i][j];
                }
            }

            // draw the heat map
            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    g2.setColor(getHeatColor(E[i][j],min,max));
                    g2.fillRect(i,n-1-j,1,1);
                }
            }
        }

        private void drawBoard(Graphics2D g2, int w, int h) {
            // we assume 1 mm per pixel resolution
            // (although this could be wrong, if the layout manager messed up)

            g2.setColor(Color.BLACK);

            // outer double ring
            g2.drawOval(-1,1,w,h);

            // inner double ring
            g2.drawOval((int)(Stats.R-1-Stats.R5),(int)(Stats.R+1-Stats.R5),
                    (int)(2*Stats.R5),(int)(2*Stats.R5));

            // outer triple ring
            g2.drawOval((int)(Stats.R-1-Stats.R4),(int)(Stats.R+1-Stats.R4),
                    (int)(2*Stats.R4),(int)(2*Stats.R4));

            // inner triple ring
            g2.drawOval((int)(Stats.R-1-Stats.R3),(int)(Stats.R+1-Stats.R3),
                    (int)(2*Stats.R3),(int)(2*Stats.R3));

            // outer bullseye ring
            g2.drawOval((int)(Stats.R-1-Stats.R2),(int)(Stats.R+1-Stats.R2),
                    (int)(2*Stats.R2),(int)(2*Stats.R2));

            // inner bullseye ring
            g2.drawOval((int)(Stats.R-1-Stats.R1),(int)(Stats.R+1-Stats.R1),
                    (int)(2*Stats.R1),(int)(2*Stats.R1));

            // draw the lines
            for (int i=0; i<10; i++) {
                double th = Math.PI*11/20 - i*Math.PI/10;
                g2.drawLine((int)Math.round(Stats.R-1+Stats.R*Math.cos(th)),
                        (int)Math.round(Stats.R+1-Stats.R*Math.sin(th)),
                        (int)Math.round(Stats.R-1+Stats.R2*Math.cos(th)),
                        (int)Math.round(Stats.R+1-Stats.R2*Math.sin(th)));
                g2.drawLine((int)Math.round(Stats.R-1+Stats.R*Math.cos(th+Math.PI)),
                        (int)Math.round(Stats.R+1-Stats.R*Math.sin(th+Math.PI)),
                        (int)Math.round(Stats.R-1+Stats.R2*Math.cos(th+Math.PI)),
                        (int)Math.round(Stats.R+1-Stats.R2*Math.sin(th+Math.PI)));
            }

            // draw the numbers
            FontMetrics fm = g2.getFontMetrics();
            int oy = fm.getAscent()/2;
            for (int i=0; i<20; i++) {
                double th = Math.PI/2 - i*Math.PI/10;
                String str = Integer.toString(Stats.d[i]);
                int ox = fm.stringWidth(str)/2;
                g2.drawString(str,
                        (int)(Stats.R-1+(Stats.R-25)*Math.cos(th)-ox),
                        (int)(Stats.R+1-(Stats.R-25)*Math.sin(th)+oy));
            }

        }

        private void drawOptimum(Graphics2D g2, int w, int h) {
            // we assume 1 mm per pixel

            g2.setColor(Color.BLACK);
            g2.fillOval(DartsApplet.this.imax-3,h-1-DartsApplet.this.jmax-3,6,6);
        }

        public void clear() {
            bufferedImage = null;
        }
    }

    private void setColorScheme(Color[] colorScheme) {
        heatMap.setColorScheme(colorScheme);
        legendLabel.setIcon(new ImageIcon(heatMap.makeLegend(100,30)));
    }

    private static final int INVALID_FORMAT_FAIL = 1;
    private static final int INVALID_NUMBERS_FAIL = 2;
    private static final int GENERAL_FAIL = 3;
    private static final String[] failMessages = new String[] {"",
            "There is a problem with the format of your scores.",
            "The following scores are not valid: ",
            "Something went wrong when computing the heat map."};
    private static String invalidNumbers;

    private StatsTask task;

    private int createHeatMap() {
        String str = textArea.getText();
        if (str == null) return INVALID_FORMAT_FAIL;

        // try to parse the text to get integer scores
        String[] scoresStr = str.split(",");
        int n = scoresStr.length;
        int[] scores = new int[n];
        try {
            for (int i=0; i<n; i++) {
                scores[i] = Integer.parseInt(scoresStr[i].trim());
            }
        }
        catch (NumberFormatException e) {
            // we failed
            return INVALID_FORMAT_FAIL;
        }
        if (scores == null) return INVALID_FORMAT_FAIL;

        // check for invalid numbers (unachievable scores)
        checkNumbers(scores);
        if (invalidNumbers.length() > 0) return INVALID_NUMBERS_FAIL;

        // everything passed, so make a new task
        task = new StatsTask(scores);
        task.start();

        return 0;
    }

    private static void checkNumbers(int[] scores) {
        Set<Integer> invalids = new TreeSet<Integer>();
        boolean atCapacity = false;

        for (int i=0; i<scores.length; i++) {
            int s = scores[i];
            if (((0 <= s) && (s <= 20)) || s == 21 || s == 22 || s == 24 || s == 25 || s == 26 ||
                    s == 27 || s == 28 || s == 30 || s == 32 || s == 33 || s == 34 || s == 36 ||
                    s == 38 || s == 39 || s == 40 || s == 42 || s == 45 || s == 48 || s == 50 ||
                    s == 51 || s == 54 || s == 57 || s == 60) {
                continue;
            }
            else {
                if (invalids.size() == 10) {
                    atCapacity = true;
                    break;
                }
                else invalids.add(s);
            }
        }

        Iterator<Integer> it = invalids.iterator();
        StringBuffer sb = new StringBuffer();
        while (it.hasNext()) {
            sb.append(Integer.toString(it.next())+", ");
        }

        if (sb.length() > 0) {
            // remove the last ", "
            sb.delete(sb.length()-2,sb.length());

            if (atCapacity) sb.append(" ...");
            else sb.append(".");
        }

        invalidNumbers = sb.toString();
    }

    private double[][] E;
    private double s1, s2, s3;
    private int imax, jmax;

    private class StatsTask extends Thread {
        private int[] scores;

        public StatsTask(int[] scores) {
            this.scores = new int[scores.length];
            for (int i=0; i<scores.length; i++) {
                this.scores[i] = scores[i];
            }
        }

        public void run() {
            enableAll(false);

            // our resolution: 1 mm per pixel)
            int nn = (int)(2*Stats.R);
            E = null;

            if (simpleModel) {
                // we're going to use the simple model

                // run the EM algorithm
                updateProgress("<html>Running EM<br>algorithm...</html>");
                s1 = Stats.simpleEM(scores,DartsApplet.this);

                // compute the heat map
                updateProgress("<html>Computing expected<br>scores...</html>");
                E = Stats.computeExpScores(s1,nn,DartsApplet.this);
            }
            else {
                // we're going to use the general model

                // run the EM algorithm
                updateProgress("<html>Running EM<br>algorithm...</html>");
                double[] a = Stats.generalEM(scores,DartsApplet.this);
                s1=a[0]; s2=a[1]; s3=a[2];

                // compute the heat map
                updateProgress("<html>Computing expected<br>scores...</html>");
                E = Stats.computeExpScores(s1,s2,s3,nn,DartsApplet.this);
            }

            // tell the heat map to draw itself
            updateProgress("<html>Drawing the<br>heat map...</html>");
            heatMap.clear();
            heatMap.repaint();

            // update all the labels
            double[] b = Stats.getMaxAndArgmax(E);
            double max = b[0];
            imax = (int)b[1]; jmax = (int)b[2];
            String aimStr = Stats.getRegion(imax,jmax,nn,nn);
            int n = scores.length;

            maxExpScoreField.setText(Double.toString(roundTo2(max)));
            aimForField.setText(aimStr);
            if (simpleModel) {
                sigmaXField.setText(Double.toString(roundTo2(Math.sqrt(s1))));
                sigmaYField.setText(" ");
                rhoField.setText(" ");
            }
            else {
                sigmaXField.setText(Double.toString(roundTo2(Math.sqrt(s1))));
                sigmaYField.setText(Double.toString(roundTo2(Math.sqrt(s2))));
                rhoField.setText(Double.toString(roundTo2(s3/Math.sqrt(s1*s2))));
            }
            numScoresLabel.setText("<html>Computed based on<br>" +
                    n + " score" + ((n>1)?"s":"") + ".</html>");

            updateProgress(112);
            updateProgress("<html>&nbsp;<br>Done.</html>");

            enableAll(true);
        }
    }

    private void enableAll(boolean b) {
        createHeatMapButton.setEnabled(b);
        changeColorsButton.setEnabled(b);
        showBoardCheckBox.setEnabled(b);
        showArgmaxCheckBox.setEnabled(b);
        simpleButton.setEnabled(b);
        generalButton.setEnabled(b);
        if (b) heatMap.addMouseMotionListener(this);
        if (!b) heatMap.removeMouseMotionListener(this);
    }

    public void updateProgress(int progress) {
        progressBar.setValue(progress);
    }

    public void updateProgress(String progressString) {
        progressLabel.setText(progressString);
    }

    public void start() {
    }

    public void stop() {
    }

    public void destroy() {
    }

    public static void main(String[] args) {
        Runnable guiCreator = new Runnable() {
            public void run() {
                JFrame window = new JFrame("A Statistician Plays Darts");
                window.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

                new DartsApplet().init(window);

                window.setSize(800, 600);
                window.setVisible(true);
            }
        };

        SwingUtilities.invokeLater(guiCreator);
    }
}
