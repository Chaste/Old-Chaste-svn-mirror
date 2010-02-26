/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.StringTokenizer;
import java.lang.Math;
import javax.imageio.ImageIO;

import javax.swing.JPanel;
import javax.swing.JLabel;

import java.text.DecimalFormat;
import java.text.NumberFormat;

public class Visualize2dCells implements ActionListener, AdjustmentListener, ItemListener, Runnable
{
    private Thread updateThread;

    static CustomCanvas2D canvas;

    Button run;
    
    public Frame frame = new Frame();
    
    public static boolean parsed_all_files = false;
    public static boolean drawAncestors = false;
    public static boolean drawAxes = true;
    public static boolean drawCells = true;
    public static boolean drawSprings = false;
    public static boolean drawCircles = false;
    public static boolean drawNutrient = false;
    public static boolean drawBetaCatenin = false;
    public static boolean drawAverageStress = false;
    public static boolean drawDifferenceStress = false;
    public static boolean writeFiles = false;
    public static boolean drawGhosts = false;
    public static boolean drawFibres = false;
    public static boolean drawCylinder = false;
    public static boolean drawCylinderOverride = true;
    public static boolean ancestorsFilePresent = false;
    public static boolean fibresFilePresent = false;
    public static boolean nutrientFilePresent = false;
    public static boolean orientationFilePresent = false;
    public static boolean betaCateninFilePresent = false;
    public static boolean stressFilePresent = false;
    public static boolean elementFilePresent = true;
    public static boolean isSparseMesh = false;
    public static boolean drawOrientations = false;
    // by default the last timestep isn't read or visualised; this
    // allows the visualiser to be run as a simulation is being run.
    // To visualise the last timestep, use "showlaststep" as an argument
    public static boolean showLastStep = false; 
    public static boolean axesEqual = false;
    
    public static int timeStep = 0;
    public static int delay = 50;
    public static int numSteps = 0;
    public static int memory_factor = 2;
    public static int[] numCells;
    public static int[] numElements;
    public static int[][] ancestor_values;
    public static int[][] element_nodes;
    public static int[][] cell_type;
    public static int[][] image_cells;
        
    public static double max_x = -1e12;
    public static double max_y = -1e12;
    public static double min_x =  1e12;
    public static double min_y =  1e12;
    public static double crypt_width = 0.0;
    public static double half_width = 0.0;
    public static double stress_time = 0.0;
    public static double force_cutoff = 0.0;
    public static double[] times;
    public static double[][] nutrient_values;
    public static double[][][] beta_catenin_values;
    public static double[][][] stress_values;
    
    public static RealPoint[][] positions;
    public static RealPoint[][] fibres;
    public static RealPoint[][] orientations;
    
    public static Scrollbar delay_slider = new Scrollbar(Scrollbar.HORIZONTAL, delay, 1, 1, 100);
    public static Scrollbar time_slider = new Scrollbar(Scrollbar.HORIZONTAL, timeStep, 1, 0, 2);
    
    public static Checkbox output = new Checkbox("Output");
    public static Checkbox springs = new Checkbox("Springs");
    public static Checkbox fibre = new Checkbox("Fibres");
    public static Checkbox cells = new Checkbox("Cells");
    public static Checkbox ghost_nodes = new Checkbox("Ghosts");
    public static Checkbox circles = new Checkbox("Cells as circles");
    public static Checkbox nutrient = new Checkbox("Nutrient");
    public static Checkbox beta_catenin = new Checkbox("Beta catenin");
    public static Checkbox average_stress = new Checkbox("Average Stress");
    public static Checkbox difference_stress = new Checkbox("Difference Stress");
    public static Checkbox ancestors = new Checkbox("Clonal Populations");
    public static Checkbox axes = new Checkbox("Axes");
    public static Checkbox axes_equal = new Checkbox("Axes Equal");
    public static Checkbox orientation = new Checkbox("Orientation");
    
    public static JLabel nearest_label = new JLabel();
    
    public static File node_file;
    public static File cell_type_file;
    public static File element_file;
    public static File nutrient_file;
    public static File beta_catenin_file;
    public static File stress_file;
    public static File ancestors_file;
    public static File fibre_file;
    public static File setup_file;
    public static File orientation_file;
    
    public static Button refresh;
    
    public Visualize2dCells() 
    {
        frame.setSize(700, 700);
        frame.setLayout(new BorderLayout());
        
        canvas = new CustomCanvas2D(this);
        canvas.setPreferredSize(new Dimension(frame.getWidth(),frame.getHeight()));
        canvas.addMouseMotionListener(canvas);
        
        addButtons(frame);
        addTimeSlider(frame);
        
        JPanel canvasPanel = new JPanel();
        canvasPanel.add(canvas);
        frame.add(canvasPanel, BorderLayout.CENTER);
                
        frame.addWindowListener(new WindowAdapter() 
        {
            public void windowClosing(WindowEvent e) 
            {
                System.exit(0);
            }
        });
        frame.pack();
        frame.setVisible(true);
    }

    public void actionPerformed(ActionEvent event) 
    {
        String pressed = event.getActionCommand();
                
        if (pressed == "Quit") 
        {
            frame.dispose();
        }
        if (pressed == "Run") 
        {
            if (timeStep == numSteps - 1) 
            {
                timeStep = 0;
                time_slider.setValue(timeStep);
            }
            if (updateThread == null) 
            {
                run.setLabel("Pause");
                updateThread = new Thread(this);
                updateThread.start();
            }
        }
        if (pressed == "Reset") 
        {
            timeStep = 0;
            time_slider.setValue(timeStep);
            canvas.drawBufferedImage();
            canvas.repaint();
        }
        if (pressed == "Pause") 
        {
            if (updateThread != null) 
            {
                Thread fred = updateThread;
                updateThread = null;
                fred.interrupt();
                run.setLabel("Run");
            }
        }
        if (pressed == "Refresh")
        {
            refresh.setEnabled(false);
        	LoadAllFiles();
            canvas.drawBufferedImage();
            canvas.repaint();
            refresh.setEnabled(true);
        }
    }
    
    public void itemStateChanged(ItemEvent e) 
    {
        Object cb = e.getItemSelectable();
        boolean state = (e.getStateChange() == ItemEvent.SELECTED);
               
        if (cb == output) 
        {
            writeFiles = state;
            System.out.println("Writing output files = "+writeFiles);
        } 
        else if (cb == springs) 
        {
            drawSprings = state;
            System.out.println("Drawing springs = "+drawSprings);
        }
        else if (cb == fibre)
        {
            drawFibres = state;
            System.out.println("Drawing fibres = "+drawFibres);
        }
        else if (cb == cells)
        {
            drawCells = state;
            if (state)
            {
            	drawBetaCatenin = false;
            	beta_catenin.setState(false);
            	drawCircles = false;
                circles.setState(false);
            }
            System.out.println("Drawing cells = "+drawCells);
        }
        else if (cb == ghost_nodes)
        {
            drawGhosts = state;
            System.out.println("Drawing ghost nodes = "+drawGhosts);    
        }
        else if (cb == circles)
        {
            drawCircles = state;
            if (state)
            {
            	drawBetaCatenin = false;
            	beta_catenin.setState(false);
            	drawCells = false;
                cells.setState(false);
            }
            System.out.println("Drawing cells as circles = "+drawCircles);    
        }
        else if (cb == nutrient)
        {
            drawNutrient = state;
            if (state)
            {
            	drawCells = false;
            	cells.setState(false);
            	drawAverageStress = false;
            	average_stress.setState(false);
            	drawDifferenceStress = false;
            	difference_stress.setState(false);
            	drawBetaCatenin = false;
            	beta_catenin.setState(false);
            }
            System.out.println("Drawing nutrient = "+drawNutrient); 
        }
        else if (cb == beta_catenin)
        {
            drawBetaCatenin = state;
            if (state)
            {
            	drawCells = false;
            	cells.setState(false);
            	circles.setState(false);
            	circles.setVisible(false);
            	drawCircles = false;
            }
            else
            {
            	circles.setVisible(true);
            }
            System.out.println("Drawing beta catenin = "+drawBetaCatenin); 
        }
        else if (cb == average_stress)
        {
            drawAverageStress = state;
            if (state)
            {
            	drawDifferenceStress = false;
            	difference_stress.setState(false);
            }
            System.out.println("Drawing average stress = "+drawAverageStress); 
        }
        else if (cb == difference_stress)
        {
            drawDifferenceStress = state;
            if (state)
            {
            	drawAverageStress = false;
            	average_stress.setState(false);
            }
            System.out.println("Drawing difference stress = "+drawDifferenceStress); 
        }
        else if (cb == axes)
        {
            drawAxes = state;
            System.out.println("Drawing axes = "+drawAxes); 
        }
        else if (cb == axes_equal)
        {
            axesEqual = state;
            
            if (axesEqual==true) // Make the axes equal
            {
	            if (min_y < min_x)
	        	{
	        		min_x = min_y;
	        	}
	            else
	            {
	            	min_y = min_x;
	            }
	            if (max_x < max_y)
	        	{
	        		max_x = max_y;
	        	}
	            else
	            {
	            	max_y = max_x;
	            }
            }
            else // reset the axes
            { 
            	CalculateCanvasDimensions();
            }
            System.out.println("Drawing axes equal = "+axesEqual); 
        }
        else if (cb == ancestors)
        {
            drawAncestors = state;
            if (state)
            {
            	drawBetaCatenin = false;
            	beta_catenin.setState(false);
            	drawAverageStress = false;
            	average_stress.setState(false);
            	drawDifferenceStress = false;
            	difference_stress.setState(false);
            }
            System.out.println("Drawing clonal populations = " + drawAncestors); 
        }
        else if (cb == orientation)
        {
            drawOrientations = state;
            System.out.println("Drawing cell orientations = " + drawOrientations); 
        }
        canvas.drawBufferedImage();
        canvas.repaint();
    }

    public void adjustmentValueChanged(AdjustmentEvent e) 
    {
        delay = delay_slider.getValue();
        timeStep = time_slider.getValue();
        canvas.drawBufferedImage();
        canvas.repaint();
    }

    public void run() 
    {
        while (updateThread != null) 
        {            
            if (timeStep < numSteps - 1) 
            {
            	timeStep++;
            } 
            else 
            {
                if (updateThread != null) 
                {
                    Thread thread = updateThread;
                    updateThread = null;
                    thread.interrupt();
                    run.setLabel("Run");
                }
            }
            try 
            {
                Thread.sleep((100 - delay) * 1);
            } 
            catch (InterruptedException e) 
            {
                return;
            }
            canvas.drawBufferedImage();
            canvas.repaint();
        }
    }

    public void addButtons(Frame frame) 
    {
        JPanel buttonPanel = new JPanel(new GridLayout(0,4));
        Button quit = new Button("Quit");
        quit.addActionListener(this);

        run = new Button("Run");
        run.addActionListener(this);
        
        Button reset = new Button("Reset");
        reset.addActionListener(this);
        
        refresh = new Button("Refresh");
        refresh.setEnabled(true);
        refresh.addActionListener(this);
        
        buttonPanel.add(quit);
        buttonPanel.add(run);
        buttonPanel.add(reset);
        buttonPanel.add(refresh);
                
        JPanel scrollPanel = new JPanel();
        
        delay_slider.setPreferredSize(new Dimension(frame.getWidth(),20));
        delay_slider.addAdjustmentListener(this);
        
        Label slow = new Label("Slow");
        Label fast = new Label("Fast");
        
        scrollPanel.add(slow);
        scrollPanel.add(delay_slider);
        scrollPanel.add(fast);
                
        JPanel northPanel = new JPanel(new GridLayout(2,0));
        northPanel.add(buttonPanel);
        northPanel.add(scrollPanel);
        
        frame.add(northPanel,BorderLayout.NORTH);
    }
    
    public void addTimeSlider(Frame frame) 
    {
        JPanel scrollPanel_time = new JPanel();
        
        time_slider.setPreferredSize(new Dimension(frame.getWidth(),20));
        time_slider.addAdjustmentListener(this);
        
        Label start_time = new Label("Start");
        Label end_time = new Label("End");
        
        scrollPanel_time.add(start_time);
        scrollPanel_time.add(time_slider);
        scrollPanel_time.add(end_time);
                        
        JPanel checkPanel = new JPanel(new GridLayout(0,5));
        output.addItemListener(this);
        springs.addItemListener(this);
        fibre.addItemListener(this);
        cells.addItemListener(this);
        ghost_nodes.addItemListener(this);
        circles.addItemListener(this);
        nutrient.addItemListener(this);
        beta_catenin.addItemListener(this);
        average_stress.addItemListener(this);
        difference_stress.addItemListener(this);
        axes.addItemListener(this);
        axes_equal.addItemListener(this);
        ancestors.addItemListener(this);
        orientation.addItemListener(this);
        
        
        checkPanel.add(output);
        checkPanel.add(springs);
        checkPanel.add(cells);
        checkPanel.add(ghost_nodes);
        checkPanel.add(circles);
        checkPanel.add(axes);
        checkPanel.add(axes_equal);
        checkPanel.add(orientation);
        checkPanel.add(nutrient);
        checkPanel.add(beta_catenin);
        checkPanel.add(fibre);
        checkPanel.add(average_stress);
        checkPanel.add(difference_stress);
        checkPanel.add(ancestors);
        checkPanel.add(nearest_label);
        
        
        JPanel southPanel = new JPanel(new GridLayout(2,0));
        southPanel.add(scrollPanel_time);
        southPanel.add(checkPanel);
        frame.add(southPanel,BorderLayout.SOUTH);
    }    
    
    public static void main(String args[]) 
    {
        System.out.println("Copyright The Chaste Project");

        // Set default states for options
        cells.setState(true);
        output.setState(false);
        springs.setState(false);
        fibre.setState(false);
        ghost_nodes.setState(false);
        circles.setState(false);
        nutrient.setState(false);
        beta_catenin.setState(false);
        average_stress.setState(false);
        difference_stress.setState(false);
        axes.setState(true);
        axes_equal.setState(false);
        ancestors.setState(false);
        orientation.setState(false);
        
        // Update states for options according to input args
        for (int i=1; i<args.length; i++)
        {
            if (args[i].equals("output"))
            {
                writeFiles = true;
                output.setState(true);                  
            } 
            else if (args[i].equals("springs"))
            {
                drawSprings = true;
                springs.setState(true);             
            }   
            else if (args[i].equals("fibres"))
            {
                drawFibres = true;
                fibre.setState(true);               
            }   
            else if (args[i].equals("nocells"))
            {
                drawCells = false;
                cells.setState(false);
            }   
            else if (args[i].equals("ghosts"))
            {
                drawGhosts = true;
                ghost_nodes.setState(true);
            }
            else if (args[i].equals("circles"))
            {
                drawCircles = true;
                circles.setState(true);
            }
            else if (args[i].equals("nutrient"))
            {
                drawNutrient = true;
                nutrient.setState(true);
                drawCells = false;
            	cells.setState(false);
            }
            else if (args[i].equals("betacatenin"))
            {
            	drawBetaCatenin = true;
            	drawCells = false;
            	cells.setState(false);
                beta_catenin.setState(true);
            }
            else if (args[i].equals("notcylindrical"))
            {
                drawCylinderOverride = false;
            }
            else if (args[i].equals("axes"))
            {
                drawAxes = true;
                axes.setState(true);
            }
            else if (args[i].equals("axesequal"))
            {
                axesEqual = true;
                axes_equal.setState(true);
            }
            else if (args[i].equals("ancestors"))
            {
                drawAncestors = true;
                ancestors.setState(true);
            }
            else if (args[i].equals("orientation"))
            {
                drawOrientations = true;
                orientation.setState(true);
            }
            else
            {
                System.out.println("Input option not recognised");
            }
        }

        // Read in results files for visualization

        node_file = new File(args[0]+"/results.viznodes");
        cell_type_file = new File(args[0]+"/results.vizcelltypes");
        element_file = new File(args[0]+"/results.vizelements");
        nutrient_file = new File(args[0]+"/results.viznutrient");
        beta_catenin_file = new File(args[0]+"/results.vizbetacatenin");
        stress_file = new File(args[0]+"/results.vizstress");
        ancestors_file = new File(args[0]+"/results.vizancestors");
        fibre_file = new File(args[0]+"/results.vizfibres");
        orientation_file = new File(args[0]+"/results.vizorientation");
                
        if (!node_file.isFile())
        {
            System.out.println("The file "+args[0]+"/results.viznodes doesn't exist");
            return;
        }
        if (!cell_type_file.isFile())
        {
            System.out.println("The file "+args[0]+"/results.vizcelltypes doesn't exist");
            return;
        }
        if (!element_file.isFile())
        {
        	// If the results.vizelements does not exist, then assume the results 
        	// were generated using a NodeBasedTissue
        	elementFilePresent = false;
        }
        if (!nutrient_file.isFile())
        {
            nutrient.setVisible(false);
            nutrient.setState(false);
            drawNutrient = false;
        }
        else
        {
        	nutrientFilePresent = true;
        }
        if (!beta_catenin_file.isFile())
        {
            beta_catenin.setVisible(false);
            beta_catenin.setState(false);
            drawBetaCatenin = false;
        }
        else
        {
        	betaCateninFilePresent = true;
        }
        if (!stress_file.isFile())
        {
            average_stress.setVisible(false);
            difference_stress.setVisible(false);
            average_stress.setState(false);
            difference_stress.setState(false);
            drawAverageStress = false;
            drawDifferenceStress = false;
        }
        else
        {
            average_stress.setState(true);
            difference_stress.setState(false);
            drawAverageStress = true;
            drawDifferenceStress = false;
            stressFilePresent = false;
        }
        
        if (!ancestors_file.isFile())
        {
            ancestors.setVisible(false);
            ancestors.setState(false);
            drawAncestors = false;
        }
        else
        {
            ancestorsFilePresent = true;
            ancestors.setState(true);
            drawAncestors = true;
        }
       
        if (!fibre_file.isFile())
        {
            fibre.setVisible(false);
            drawFibres = false;
        } 
        else 
        {
            fibre.setState(true);
            drawFibres = true;
            fibresFilePresent = true;
        }
        
        if (!orientation_file.isFile())
        {
            orientation.setVisible(false);
            orientation.setState(false);
            drawOrientations = false;
        }
        else
        {
        	orientationFilePresent = true;
        	orientation.setState(true);
            drawOrientations = true;
        }
        
        
        setup_file = new File(args[0]+"/results.vizsetup");
        if (!setup_file.isFile())
        {
            System.out.println("The file "+args[0]+"/results.vizsetup doesn't exist");
        }
        else 
        {
            try
            {
                BufferedReader in_setup_file = new BufferedReader(new FileReader(setup_file));
                String line_setup = in_setup_file.readLine();  
                
                // Read setup information
                while (line_setup != null)
                {
                    StringTokenizer st_setup = new StringTokenizer(line_setup);
                    String parameter = st_setup.nextToken();
                    if (parameter.equals("MeshWidth"))  
                    {
                        crypt_width = Double.valueOf(st_setup.nextToken()).doubleValue();
                        half_width = crypt_width/2.0;
                        System.out.println("Mesh Width = " + crypt_width);
                        if (crypt_width < 1.0)
                        {
                        	System.out.println("Mesh Width less than 1 so can't be periodic - setting mesh to be noncylindrical");	
                            drawCylinderOverride = false;
                        }
                        drawCylinder = true && drawCylinderOverride;    // this is made true only if mesh width exists
                    }
                    if (parameter.equals("Nutrient"))  
                    {
                        // Overrule the previous bit since for nutrient sims, we don't want cylindrical periodicity
                        drawCylinder = false; 
                        drawNutrient = true;
                        nutrient.setState(true);
                    }
                    if (parameter.equals("BetaCatenin")) 
                    {
                        drawBetaCatenin = true;
                        beta_catenin.setState(true);
                        drawCells = false;
                        cells.setState(false);
                        circles.setVisible(false);
                    }
                    if (parameter.equals("Cutoff"))
                    {
                    	force_cutoff = Double.valueOf(st_setup.nextToken()).doubleValue();
                    	System.out.println("Force cutoff = " + force_cutoff);
                    }
                    if (parameter.equals("Complete")) 
                    {
                	    showLastStep = true;
                    }
                    if (parameter.equals("SparseMesh")) 
                    {
                	    isSparseMesh = true;
                    }
                    line_setup = in_setup_file.readLine();
                }
            }
            catch (Exception e) 
        	{
            	System.out.println("Error occured. Exception message: "+e.getMessage());
        	}
        }             
        
        final Visualize2dCells vis = new Visualize2dCells();

        LoadAllFiles();
        
        canvas.drawBufferedImage();
        canvas.repaint();
    }

    
    public static void LoadAllFiles()
    {
    	parsed_all_files = false;
        try 
        {                  	
        	BufferedReader skim_node_file = new BufferedReader(new FileReader(node_file));

            int num_lines = 0;
            while (skim_node_file.readLine() != null) 
            {
                num_lines++;
            }
            
            // By default we don't read or print the final line, so
            // the visualiser can be run as a simulation is being run,
            // and the visualiser will work on incomplete data.
            if (num_lines>1 && !showLastStep)
            {
            	num_lines -= 1;
            }

            numSteps = num_lines;
            time_slider.setMaximum(numSteps);
            times = new double[num_lines];
            positions = new RealPoint[num_lines][];
            cell_type = new int [num_lines][];
            numCells = new int[num_lines];            
            image_cells = new int[num_lines][];
            
            if (elementFilePresent)
            {
            	numElements = new int[num_lines];
            	element_nodes = new int[num_lines][];
            }
            if (fibresFilePresent)
            {
            	fibres = new RealPoint[num_lines][];
            }
            if (nutrientFilePresent)
            {
            	nutrient_values = new double[num_lines][];
            }
            if (betaCateninFilePresent)
            {
            	beta_catenin_values = new double[num_lines][][];
            }
            if (stressFilePresent)
            {
            	stress_values = new double[num_lines][][];
            }
            if (ancestorsFilePresent)
            {
            	ancestor_values = new int[num_lines][];
            }
            if (orientationFilePresent)
            {
            	orientations = new RealPoint[num_lines][];
            }
            
            String line_fibre = "";
            BufferedReader in_fibre_file = null;
            if (drawFibres)
            {
                fibres = new RealPoint[num_lines][]; 
                in_fibre_file = new BufferedReader(new FileReader(fibre_file));
                line_fibre = in_fibre_file.readLine();
            }

            String line_orientation = "";
            BufferedReader in_orientation_file = null;
            if (drawOrientations)
            {
                orientations = new RealPoint[num_lines][]; 
                in_orientation_file = new BufferedReader(new FileReader(orientation_file));
                line_orientation = in_orientation_file.readLine();
            }
            
            String line_nutrient = "";
            BufferedReader in_nutrient_file = null;
            if (drawNutrient)
            {
               	nutrient_values = new double[num_lines][]; 
               	in_nutrient_file = new BufferedReader(new FileReader(nutrient_file));
               	line_nutrient = in_nutrient_file.readLine();
            }
           	
            String line_ancestors = "";
            BufferedReader in_ancestors_file = null;
            if (drawAncestors)
            {
               	ancestor_values = new int[num_lines][]; 
               	in_ancestors_file = new BufferedReader(new FileReader(ancestors_file));
               	line_ancestors = in_ancestors_file.readLine();
            }

            String line_beta_catenin = "";
            BufferedReader in_beta_catenin_file = null;
            if (drawBetaCatenin)
            {
               	beta_catenin_values = new double[num_lines][][]; 
               	in_beta_catenin_file = new BufferedReader(new FileReader(beta_catenin_file));
               	line_beta_catenin = in_beta_catenin_file.readLine();
            }
          
            String line_stress = "";
            BufferedReader in_stress_file = null;
            if (drawAverageStress || drawDifferenceStress)
            {
               	stress_values = new double[num_lines][][]; 
               	in_stress_file = new BufferedReader(new FileReader(stress_file));
               	line_stress = in_stress_file.readLine();
            }
            
            BufferedReader in_node_file = new BufferedReader(new FileReader(node_file));
            String line_node = in_node_file.readLine(); // from console input example
            
            BufferedReader in_cell_type_file = new BufferedReader(new FileReader(cell_type_file));
            String line_cell_type = in_cell_type_file.readLine(); // from console input example

            BufferedReader in_element_file = null;
            String line_element = null;
            if (elementFilePresent)
            {
            	in_element_file = new BufferedReader(new FileReader(element_file));
            	line_element = in_element_file.readLine();   // above
            }
            else
            {
            	System.out.println("No vizelement file found - plotting node locations only");
            }
            
            // If line is not end of file continue
            boolean has_stress_line_been_read = false;
            int row = 0;
            while (line_node != null && row<num_lines) 
            {
            	// Create a StringTokenizer with a colon sign as a delimiter
                StringTokenizer st_node = new StringTokenizer(line_node);
                StringTokenizer st_cell_type = new StringTokenizer(line_cell_type);
                StringTokenizer st_element = null;
                if (elementFilePresent)
                {
                	st_element = new StringTokenizer(line_element);
                }
                StringTokenizer st_fibre = null;
                StringTokenizer st_nutrient = null;
                StringTokenizer st_ancestors = null;
                StringTokenizer st_beta_catenin = null;
                StringTokenizer st_stress = null;
                StringTokenizer st_orientation = null;
                
                
                Double time = Double.valueOf(st_node.nextToken());
                Double cell_type_time = Double.valueOf(st_cell_type.nextToken());
                
                if (drawFibres)
                {
                    st_fibre = new StringTokenizer(line_fibre);
                    Double fibre_time = Double.valueOf(st_fibre.nextToken());
                    //fibre_time unused
                }
                
                
                if (drawOrientations)
                {
                    st_orientation = new StringTokenizer(line_orientation);
                    Double orientation_time = Double.valueOf(st_orientation.nextToken());
                    //orientation_time unused
                }
                if (drawNutrient)
                {
                    st_nutrient = new StringTokenizer(line_nutrient);
                    Double nutrient_time = Double.valueOf(st_nutrient.nextToken());
                    
                    int nutrient_entries = st_nutrient.countTokens();
                    if (nutrient_entries%4 != 0)
                    {
                    	System.out.println("Warning: Results from time "+time.doubleValue()+" will not be plotted as the corresponding line of the nutrient file is not of the required form: time,index,x,y,nutrient,index,x,y,nutrient,...");
                    	break;
                    }
                }
                
                if (drawAncestors)
                {
                    st_ancestors = new StringTokenizer(line_ancestors);
                    Double ancestors_time = Double.valueOf(st_ancestors.nextToken());
                }
                
                if (drawBetaCatenin)
                {
                    st_beta_catenin = new StringTokenizer(line_beta_catenin);
                    Double beta_catenin_time = Double.valueOf(st_beta_catenin.nextToken());
                    
                    // Count the number of entries in the bcat file to get num non ghosts and check correct 
                    int beta_catenin_entries = st_beta_catenin.countTokens();
	            
                    if (beta_catenin_entries%6 != 0)
                    {
                    	System.out.println("Warning: Results from time "+time.doubleValue()+" will not be plotted as the corresponding line of the beta catenin file is not of the required form: time,index,x,y,bCat_mem,bCat_cyto,bCat_nuc,index,x,y,bCat_mem,bCat_cyto,bCat_nuc,...");
                    	break;
                    }
                }
                
                if ((drawAverageStress || drawDifferenceStress) && !has_stress_line_been_read)
                {
                    st_stress = new StringTokenizer(line_stress);
                    stress_time = Double.valueOf(st_stress.nextToken()).doubleValue();
                    
                    // Count the number of entries in the bcat file to get num non ghosts and check correct 
		            int stress_entries = st_stress.countTokens();
		            
		            if (stress_entries%5 != 0)
		            {
		                throw new Exception("The stress file must take the form: time,index,x,y,min_stress,max_stress,index,x,y,min_stress,max_stress,...");
		            }
		            System.out.println("The stress file is for time = " + stress_time);
		            has_stress_line_been_read = true;
                }
                
                if (elementFilePresent)
                {
                    double element_time = Double.valueOf(st_element.nextToken()).doubleValue();
                    
                    if (Math.abs(time.doubleValue() - element_time) > 1e-6) 
                    {
                    	throw new Exception("Error: The time corresponding to each line of the element file must match that of the node file");
                    }               	
                }
                
                times[row] = time.doubleValue();

                // Count the number of entries in the node file and check correct 
                int entries = st_node.countTokens();
                                
                if (entries%2 != 0)
                {
                	System.out.println("Warning: Results from time "+time.doubleValue()+" will not be plotted as the corresponding line of the node file is not of the required form: time,x,y,x,y,...");
                	break;
                }
                numCells[row] = entries/2; 

				int entries_from_cell_type_file = st_cell_type.countTokens();
				if (numCells[row] != entries_from_cell_type_file)
				{
					System.out.println("Warning: At time "+time.doubleValue()+", node file gives "+numCells[row]+" cells, but cell type file gives "+entries_from_cell_type_file+" cells");
					break;
				}

                if (elementFilePresent)
                {
                	
                    // Count the number of entries in the element file and check correct 
                    entries = st_element.countTokens();
                    
                    if (isSparseMesh)
                    {
                    	if (entries%2 != 0)
                        {
                        	System.out.println("Warning: Results from time "+time.doubleValue()+" will not be plotted as the corresponding line of the element file is not of the required form: time,n1,n2,n1,n2..");
                        	break;
                        }
                        numElements[row] = entries/2;
                        element_nodes[row] = new int[memory_factor*2*numElements[row]];
                    }
                    else
                    {	
	                    if (entries%3 != 0)
	                    {
	                    	System.out.println("Warning: Results from time "+time.doubleValue()+" will not be plotted as the corresponding line of the element file is not of the required form: time,n1,n2,n3,n1,n2,n3..");
	                    	break;
	                    }
	                    numElements[row] = entries/3;
	                    element_nodes[row] = new int[memory_factor*3*numElements[row]];
                    }
                }                

                positions[row] = new RealPoint[memory_factor*numCells[row]];
                if (fibresFilePresent)
                {
                	fibres[row] = new RealPoint[memory_factor*numCells[row]];
                }
                
                if (orientationFilePresent)
                {
                	orientations[row] = new RealPoint[memory_factor*numCells[row]];
                }
                
                if (nutrientFilePresent)
                {
                	nutrient_values[row] = new double[memory_factor*numCells[row]];
                }
                
                if (ancestorsFilePresent)
                {
                	ancestor_values[row] = new int[memory_factor*numCells[row]];
                }
                
                if (betaCateninFilePresent)
                {
                	beta_catenin_values[row] = new double[memory_factor*numCells[row]][3];
                }
                
                if (stressFilePresent)
                {
                	stress_values[row] = new double[memory_factor*numCells[row]][2];
                }
                
                cell_type[row] = new int[memory_factor*numCells[row]];
                        
                for (int i=0; i<numCells[row]; i++) 
                {
                    double d1 = Double.valueOf(st_node.nextToken()).doubleValue();
                    double d2 = Double.valueOf(st_node.nextToken()).doubleValue();

                    if (drawFibres)
                    {
                        double f1 = Double.valueOf(st_fibre.nextToken()).doubleValue();
                        double f2 = Double.valueOf(st_fibre.nextToken()).doubleValue();
                        fibres[row][i] = new RealPoint(f1,f2);
                    }
                    
                    if (drawOrientations)
                    {
                        double o1 = Double.valueOf(st_orientation.nextToken()).doubleValue();
                        double o2 = Double.valueOf(st_orientation.nextToken()).doubleValue();
                        orientations[row][i] = new RealPoint(o1,o2);
                    }
                    
                    cell_type[row][i] = Integer.parseInt(st_cell_type.nextToken());

                    if ( cell_type[row][i] < 0 )
                    {
                        System.out.println("Error: Cell type must be a non-negative integer");
                        System.exit(0);
                    }
                    positions[row][i] = new RealPoint(d1,d2);
                    
                    if (drawNutrient)
                    {
                    	if (cell_type[row][i]!=canvas.INVISIBLE_COLOUR)	// if this is not a ghost cell
                    	{
                    		String skip; // skips past unnecessary information
                        	int index = Integer.parseInt(st_nutrient.nextToken()); // index
                        	skip = st_nutrient.nextToken(); // x
                        	skip = st_nutrient.nextToken(); // y
                        	
                        	double nutrient = Double.valueOf(st_nutrient.nextToken()).doubleValue();
                        	nutrient_values[row][index] = nutrient;
                        }
                    }	
                    
                    if (drawAncestors)
                    {	
                    	// If this is a real cell then read in ancestor from row
                    	if (cell_type[row][i]!=canvas.INVISIBLE_COLOUR)	// if this is not a ghost cell
                    	{
                    		int ancestor_value = Integer.parseInt(st_ancestors.nextToken()); // index
                        	ancestor_values[row][i] = ancestor_value;
                        }
                    }	
                    
                    if (drawBetaCatenin)
                    {
                    	if (cell_type[row][i] != canvas.INVISIBLE_COLOUR) // if this is not a ghost cell
                    	{
                    		String skip; 
                        	int index = Integer.parseInt(st_beta_catenin.nextToken()); // index
                        	skip = st_beta_catenin.nextToken(); // x
                        	skip = st_beta_catenin.nextToken(); // y
                        	
                        	double beta_catenin_membrane = Double.valueOf(st_beta_catenin.nextToken()).doubleValue();
                        	double beta_catenin_cytoplasm = Double.valueOf(st_beta_catenin.nextToken()).doubleValue();
                        	double beta_catenin_nuclear = Double.valueOf(st_beta_catenin.nextToken()).doubleValue();
                        	beta_catenin_values[row][index][0] = beta_catenin_membrane;
                        	beta_catenin_values[row][index][1] = beta_catenin_cytoplasm;
                        	beta_catenin_values[row][index][2] = beta_catenin_nuclear;
                        }
                    }	
                    
                    if (Math.abs(stress_time - times[row]) < 1e-7)
                    {
                        // There is currently only ever one timestep's worth of data in 
                        // the stress file, so this is only called when it is the correct time 
                    	if (drawAverageStress || drawDifferenceStress)
                    	{
                    		String skip; 
                        	int index = Integer.parseInt(st_stress.nextToken()); // index        	                
                        	skip = st_stress.nextToken(); // x
                        	skip = st_stress.nextToken(); // y
                        	
                        	double stress_min = Double.valueOf(st_stress.nextToken()).doubleValue();
                        	double stress_max = Double.valueOf(st_stress.nextToken()).doubleValue();
                        	stress_values[0][index][0] = 0.5*(stress_min + stress_max);
                        	stress_values[0][index][1] = 0.5*(stress_max - stress_min);                    	
                    	}
                    }
                }                
                
                if (elementFilePresent)
                {
                	if (isSparseMesh)
                    {
	                	for (int i = 0; i < 2*numElements[row]; i++) 
	                    {
	                        int node = Integer.parseInt(st_element.nextToken());
	                        element_nodes[row][i] = node;
	                    }
                    }
                	else
                	{
	                	for (int i = 0; i < 3*numElements[row]; i++) 
	                    {
	                        int node = Integer.parseInt(st_element.nextToken());
	                        element_nodes[row][i] = node;
	                    }
                	}
                }
                
                // Read next line of each file

                line_node = in_node_file.readLine();
                line_cell_type = in_cell_type_file.readLine();
                
                if (drawFibres)
                {
                    line_fibre = in_fibre_file.readLine();
                }
                if (drawOrientations)
                {
                	line_orientation = in_orientation_file.readLine();
                }
                if (drawNutrient)
                {
                	line_nutrient = in_nutrient_file.readLine();
                }
                if (drawBetaCatenin)
                {
                	line_beta_catenin = in_beta_catenin_file.readLine();
                }     
                if (drawAncestors)
                {
                	line_ancestors = in_ancestors_file.readLine();
                }        
                if (elementFilePresent)
                {
                	line_element = in_element_file.readLine();
                }
                
                row++;                
            } // end while not at end of file
            
            System.out.println("Writing output files = "+writeFiles);
            System.out.println("Drawing springs = "+drawSprings);
            System.out.println("Drawing fibres = "+drawFibres);
            System.out.println("Drawing nutrient = "+drawNutrient);
            System.out.println("Drawing beta catenin = "+drawBetaCatenin);
            System.out.println("Drawing average stress = "+drawAverageStress);
            System.out.println("Drawing difference stress = "+drawDifferenceStress);
            System.out.println("Drawing cells = "+drawCells);
            System.out.println("Drawing ghost nodes = "+drawGhosts);
            System.out.println("Drawing cylindrically = "+ drawCylinder);
            System.out.println("Drawing axes = "+ drawAxes);
            System.out.println("Drawing axes equal = "+ axesEqual);
            System.out.println("Drawing clonal populations = "+ drawAncestors);
            System.out.println("Drawing orientations = "+ drawOrientations);
            System.out.println("Using sparse mesh = "+ isSparseMesh);
            
            
            if (drawCylinder) 
            {
            	ConvertCylindricalDataToPlane();
            }
            
            CalculateCanvasDimensions();
            parsed_all_files = true;
        } 
        catch (Exception e) 
        {
            System.out.println("Error occured. Exception message: "+e.getMessage());
        }
    	    	
    }
    
    
    public static void CalculateCanvasDimensions()
    {
        max_x = -1e12;
        max_y = -1e12;
        min_x =  1e12;
        min_y =  1e12;
        for (int row=0; row<numSteps; row++)
        {
            for (int i=0; i < numCells[row]; i++) 
            {
                if (positions[row][i].x > max_x) 
                {
                    max_x = positions[row][i].x;
                }   
                if (positions[row][i].y > max_y) 
                {
                    max_y = positions[row][i].y;
                } 
                if (positions[row][i].x < min_x) 
                {
                    min_x = positions[row][i].x;
                }
                if (positions[row][i].y < min_y) 
                {
                    min_y = positions[row][i].y;
                }
            }
        }
    }
    
    public static void ConvertCylindricalDataToPlane()
    {
        // Scan through each element
        for (int time_index=0; time_index<numSteps; time_index++)
        {
            image_cells[time_index] = new int[memory_factor*numCells[time_index]]; // reserve plenty of memory
            
            // Fill image_nodes with an identity map (at each time step each node maps to itself)            
            for (int i=0; i<numCells[time_index]; i++) 
            {
                image_cells[time_index][i] = i;
            }
            
            if ( (elementFilePresent) && (!isSparseMesh) )
            {
                // Draw elements first
                for (int i=0; i<numElements[time_index]; i++)
                {   
                    // What nodes are we joining up?
                    int indexA = element_nodes[time_index][3*i];
                    int indexB = element_nodes[time_index][3*i+1];
                    int indexC = element_nodes[time_index][3*i+2];
                    
                    // Find the x-co-ords of each node
                    RealPoint rA = positions[time_index][indexA];
                    RealPoint rB = positions[time_index][indexB];
                    RealPoint rC = positions[time_index][indexC];
                    
                    // Identify edges that are oversized
                    if ((Math.abs(rA.x - rB.x) > 0.75*crypt_width)
                        ||(Math.abs(rB.x - rC.x) > 0.75*crypt_width)
                        ||(Math.abs(rA.x - rC.x) > 0.75*crypt_width))
                    {
                        MakeNewImageCell(time_index,indexA);
                        MakeNewImageCell(time_index,indexB);
                        MakeNewImageCell(time_index,indexC);
                        
                        // Break those elements into two separate elements
                        SplitElement(time_index,i);
                    }
                }
            }            
        }
    }

    public static void SplitElement(int time_index,int element_index)
    {
        int indexA = element_nodes[time_index][3*element_index];
        int indexB = element_nodes[time_index][3*element_index+1];
        int indexC = element_nodes[time_index][3*element_index+2];

        // Find the x-co-ords of each node
        RealPoint rA = positions[time_index][indexA];
        RealPoint rB = positions[time_index][indexB];
        RealPoint rC = positions[time_index][indexC];
        
        // Create a new element which contains the image of node A.
        element_nodes[time_index][3*numElements[time_index]] = image_cells[time_index][indexA];
        
        // Leave Node A in this element
        
        // If node B is far away 
        if (Math.abs(rA.x - rB.x) > 0.75*crypt_width)
        {
            // Add node B to new element and add image of B to this element.
            element_nodes[time_index][3*numElements[time_index]+1] = indexB;
            element_nodes[time_index][3*element_index+1] = image_cells[time_index][indexB];
        }
        else
        {   
        	// Node B is still in this element - add its image to new element
            element_nodes[time_index][3*numElements[time_index]+1] = image_cells[time_index][indexB];
        }
        
        // If node C is far away 
        if (Math.abs(rA.x - rC.x) > 0.75*crypt_width)
        {   
        	// Add it to new element and add image of C to this element.
            element_nodes[time_index][3*numElements[time_index]+2] = indexC;
            element_nodes[time_index][3*element_index+2] = image_cells[time_index][indexC];
        }
        else
        {   
        	// Node C is still in this element - add its image to new element
            element_nodes[time_index][3*numElements[time_index]+2] = image_cells[time_index][indexC];
        }       
        numElements[time_index]++;
    }
    
    public static void MakeNewImageCell(int time_index, int node_index)
    {   
    	// Only make a new cell if one hasn't already been made
        if (image_cells[time_index][node_index] == node_index)
        {	
        	// Make a new image of Cell A       
            RealPoint new_point = positions[time_index][node_index];
            RealPoint new_point2 = new RealPoint(0.0,0.0);
            new_point2.y = new_point.y;
            
            if (new_point.x < half_width)
            {
                new_point2.x = new_point.x + crypt_width;
            }
            if (new_point.x > half_width)
            {
                new_point2.x = new_point.x - crypt_width;
            }

            // New ghost node
            positions[time_index][numCells[time_index]] = new_point2;
            cell_type[time_index][numCells[time_index]] = canvas.INVISIBLE_COLOUR;
            
            // Update the image record
            image_cells[time_index][node_index] = numCells[time_index]; 
            numCells[time_index]++;
        }    
    }
}


class RealPoint
{
    public double x, y; 
    public RealPoint(double xs, double ys)
    {
        x = xs;
        y = ys;
    }
    public RealPoint(RealPoint p1, RealPoint p2) // returns the average of p1 and p2
    {
        x = (p1.x + p2.x)/2.0;
        y = (p1.y + p2.y)/2.0;
    }
}


class PlotPoint
{
    public int x, y;
    public PlotPoint(int xs, int ys)
    {        
        x = xs;
        y = ys;
    }
}


class CustomCanvas2D extends Canvas implements MouseMotionListener 
{
    private static final long serialVersionUID = 6997195399856046957L;

    Visualize2dCells vis;

    boolean imageReady = false;   
    boolean imageDrawing = false;
    
    int width;
    int height;
    int node_radius = 2;
    
    BufferedImage buffered_image = null;
    Graphics g2 = null;
    

    Color background_white = new Color(255,255,255);
    Color spring_silver = new Color(200,200,200);
    Color apoptotic_grey = new Color(80,80,80);
    Color purple = new Color(121,126,234);
    
    public CustomCanvas2D(Visualize2dCells v) 
    {
        vis = v;
        setBackground(background_white);
    }

    public void paint(Graphics graphics)
    {
        if (vis.parsed_all_files == false)
        {
            graphics.drawString("Still parsing input...", 10,10);
            return;
        }
                
        if (!imageReady)
        {
            repaint();
        }
        imageDrawing = true;
        graphics.drawImage(buffered_image,0,0,this);
        
        if (vis.writeFiles)
        {
        	//String filename = String.format("image%1$05d.png", vis.timeStep); //Pre Java-1.6
        	String filename = String.format("image%1$05d.png", new Object[] {new Integer(vis.timeStep)});
            System.out.println("Writing file : "+filename+".");
            File f = new File(filename);
            try 
            {
                ImageIO.write(buffered_image, "png", f);
            } 
            catch (Exception e)
            {
            }
            System.out.println("Written file : "+filename+".");
        }
        imageDrawing = false;
    }
    
    public void drawBufferedImage() 
    {
        int cycle = 0;
        while (imageDrawing)
        {
            System.out.print("");
            if (cycle == 100000)
            {
                System.out.print(".");
                cycle = 0;
            }
            cycle++;
        }
        imageReady = false;
        
        int tick_length = 10;
        
        vis.time_slider.setValue(vis.timeStep); 
        
        if (g2==null)
        {
            height = getHeight();
            width = getWidth();
            buffered_image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
            g2 = buffered_image.getGraphics();
        }
        
        g2.setColor(background_white);
        g2.fillRect(0, 0, width, height);
        g2.setColor(Color.black);
        g2.drawString("Time = " + vis.times[vis.timeStep], 10, 10);

        if (vis.drawCircles)
        {
	        // Draw cell circle interiors
	        for (int i=0; i<vis.numCells[vis.timeStep]; i++ ) 
	        {
	        	SetCellColour(i); 
	            PlotPoint p = scale(vis.positions[vis.timeStep][i]);
	            int rx = (int) (0.5* width /(vis.max_x - vis.min_x));
	            int ry = (int) (0.5 * height /(vis.max_y - vis.min_y));
	            
	            if (vis.cell_type[vis.timeStep][i] != INVISIBLE_COLOUR) // if not a ghost node
	            {
	            	g2.fillOval(p.x-rx, p.y-ry, 2*rx, 2*ry);
	            }
	        }        
	        
	        // Draw cell circle boundaries
	        for (int i=0; i<vis.numCells[vis.timeStep]; i++ ) 
	        {
	        	g2.setColor(Color.black);
	            PlotPoint p=scale(vis.positions[vis.timeStep][i]);
	            int rx = (int) (0.5* width /(vis.max_x - vis.min_x));
	            int ry = (int) (0.5 * height /(vis.max_y - vis.min_y));
	            
	            if (vis.cell_type[vis.timeStep][i] != INVISIBLE_COLOUR) // if not ghost
	            {
	            	g2.drawOval(p.x-rx, p.y-ry, 2*rx, 2*ry);
	            }
	        }        
        }
        
        if (vis.drawNutrient && vis.drawCircles)
        {
            // Draw cell circle interiors
	        for (int i=0; i<vis.numCells[vis.timeStep]; i++ ) 
	        {
	        	// Don't draw ghost nodes on labelled cells
	        	if ( (vis.cell_type[vis.timeStep][i] != INVISIBLE_COLOUR) 
		            && (vis.cell_type[vis.timeStep][i] != LABELLED_COLOUR)) 	
		        {
	        		SetCellColour(i); 
	        		PlotPoint p = scale(vis.positions[vis.timeStep][i]);
	        		int rx = (int) (0.5* width /(vis.max_x - vis.min_x));
	        		int ry = (int) (0.5 * height /(vis.max_y - vis.min_y));
	        		g2.fillOval(p.x-rx, p.y-ry, 2*rx, 2*ry);
	            }
	        }
	        
	        // Draw labelled cells now, so they're on top
	        for (int i=0; i<vis.numCells[vis.timeStep]; i++ ) 
	        {
	        	if (vis.cell_type[vis.timeStep][i] == LABELLED_COLOUR)	
	        	{
	        		SetCellColour(i); 
	        		PlotPoint p = scale(vis.positions[vis.timeStep][i]);
	        		int rx = (int) (0.5* width /(vis.max_x - vis.min_x));
	        		int ry = (int) (0.5 * height /(vis.max_y - vis.min_y));
	        		g2.fillOval(p.x-rx, p.y-ry, 2*rx, 2*ry);
	        	}
	        }	        
        }

        g2.setColor(Color.black);
        Shape original_clip = g2.getClip();

        if (vis.elementFilePresent && !vis.isSparseMesh )
        {
	        // Draw elements first
	        for (int i=0; i<vis.numElements[vis.timeStep]; i++)
	        {
	            // What nodes are we joining up?
	        	int index[] = new int[3];
	            index[0] = vis.element_nodes[vis.timeStep][3*i];
	            index[1] = vis.element_nodes[vis.timeStep][3*i+1];
	            index[2] = vis.element_nodes[vis.timeStep][3*i+2];

	            RealPoint r1 = vis.positions[vis.timeStep][index[0]];
	            RealPoint r2 = vis.positions[vis.timeStep][index[1]];
	            RealPoint r3 = vis.positions[vis.timeStep][index[2]];

	            RealPoint circumcentre = DrawCircumcentre(r1, r2, r3);
	            PlotPoint plotcircumcentre = scale(circumcentre);

	            // Where are they? Convert to integer pixels
	            PlotPoint vertex[] = new PlotPoint[3];
	            vertex[0] = scale(r1);
	            vertex[1] = scale(r2);
	            vertex[2] = scale(r3);

	            PlotPoint midpoint[] = new PlotPoint[3];
	            midpoint[2] = scale(new RealPoint(r1, r2));
	            midpoint[0] = scale(new RealPoint(r2, r3));
	            midpoint[1] = scale(new RealPoint(r3, r1));

	            g2.setColor(Color.black);

	            if (vis.drawCells)
	            {
	                int clipx[] = new int[3];
	                int clipy[] = new int[3];
	                for (int node=0; node<3; node++)
	                {
	                	clipx[node] = vertex[node].x;
	                    clipy[node] = vertex[node].y;
	                }
	                Polygon clip = new Polygon(clipx, clipy, 3);
	                boolean clip_me = false;

	                // Is circumcentre in the triangle? If not, then we'll clip
	                // the next bit of drawing to fit inside the triangle
	                if (!clip.contains(new Point(plotcircumcentre.x, plotcircumcentre.y)))
	                {
	                	clip_me = true;
	                    g2.setClip(clip);
	                }
	                for (int node=0; node<3; node++)
	                {
	                    SetCellColour(index[node]);
	                    int xs[] = new int[4];
	                    int ys[] = new int[4];
	                    xs[0] = plotcircumcentre.x;
	                    ys[0] = plotcircumcentre.y;
	                    xs[1] = midpoint[(node+1)%3].x;
	                    ys[1] = midpoint[(node+1)%3].y;
	                    xs[2] = vertex[node].x;
	                    ys[2] = vertex[node].y;
	                    xs[3] = midpoint[(node+2)%3].x;
	                    ys[3] = midpoint[(node+2)%3].y;
	                    g2.fillPolygon(xs, ys, 4);
	                }

	                g2.setColor(Color.black);

	                // Plot cell boundary lines
	                if ( (vis.cell_type[vis.timeStep][index[0]]!= INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[1]]!= INVISIBLE_COLOUR) )
	                {
	                    g2.drawLine(midpoint[2].x, midpoint[2].y, plotcircumcentre.x, plotcircumcentre.y);
	                }
	                if ( (vis.cell_type[vis.timeStep][index[1]]!= INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[2]]!= INVISIBLE_COLOUR) )
	                {
	                    g2.drawLine(midpoint[0].x, midpoint[0].y, plotcircumcentre.x, plotcircumcentre.y);
	                }
	                if ( (vis.cell_type[vis.timeStep][index[2]]!= INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[0]]!= INVISIBLE_COLOUR) )
	                {
	                    g2.drawLine(midpoint[1].x, midpoint[1].y, plotcircumcentre.x, plotcircumcentre.y);
	                }
	                if (clip_me)
	                {
	                	g2.setClip(original_clip);
	                }                
	            }	        
            
	            if (vis.drawNutrient && !vis.drawCircles)
	            {                    
	            	int clipx[] = new int[3];
	            	int clipy[] = new int[3];
	            	for (int node=0; node<3; node++)
	            	{
	            		clipx[node] = vertex[node].x;
	            		clipy[node] = vertex[node].y;
	            	}
	            	Polygon clip = new Polygon(clipx, clipy, 3);
	            	boolean clip_me = false;
	            	// Is circumcentre in the triangle?
	            	// If not, then we'll clip the next bit of drawing to fit inside the triangle (ticket #432)
	            	if (!clip.contains(new Point(plotcircumcentre.x, plotcircumcentre.y)))
	            	{
	            		clip_me = true;
	            		g2.setClip(clip);
	            	}

	            	PlotPoint cutoff = scale(vis.force_cutoff, 0.0);
	            	for (int node=0; node<3; node++)
	            	{
	            		// See #542
	            		PlotPoint cutoff_point = new PlotPoint(vertex[node].x + cutoff.x, vertex[node].y + cutoff.y);
	            		int sq_dist_cutoff = SquaredDistance(vertex[node], cutoff_point);
	            		int sq_dist_prev_node = SquaredDistance(vertex[node], vertex[(node+2)%3]);
	            		int sq_dist_next_node = SquaredDistance(vertex[node], vertex[(node+1)%3]);
	                    if ( (sq_dist_prev_node <= sq_dist_cutoff) && (sq_dist_next_node <= sq_dist_cutoff))
	                    {
	                    	SetCellColour(index[node]);
	                    	int xs[] = new int[4];
	                    	int ys[] = new int[4];
	                    	xs[0] = plotcircumcentre.x;
	                    	ys[0] = plotcircumcentre.y;
	                    	xs[1] = midpoint[(node+1)%3].x;
	                    	ys[1] = midpoint[(node+1)%3].y;
	                    	xs[2] = vertex[node].x;
	                    	ys[2] = vertex[node].y;
	                    	xs[3] = midpoint[(node+2)%3].x;
	                    	ys[3] = midpoint[(node+2)%3].y;
	                    	g2.fillPolygon(xs, ys, 4);
	                    }
	            	}
	            	 	           
	            	g2.setColor(Color.black);
	            	// Plot cell boundary lines
	            	if( (vis.cell_type[vis.timeStep][index[0]] != INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[1]] != INVISIBLE_COLOUR))
	            	{
	            		g2.drawLine(midpoint[2].x, midpoint[2].y, plotcircumcentre.x, plotcircumcentre.y);
	            	}
	            	if( (vis.cell_type[vis.timeStep][index[1]] != INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[2]] != INVISIBLE_COLOUR))
	            	{
	            		g2.drawLine(midpoint[0].x, midpoint[0].y, plotcircumcentre.x, plotcircumcentre.y);
	            	}
	            	if( (vis.cell_type[vis.timeStep][index[2]] != INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[0]] != INVISIBLE_COLOUR))
	            	{
	            		g2.drawLine(midpoint[1].x, midpoint[1].y, plotcircumcentre.x, plotcircumcentre.y);
	            	}
	            	if (clip_me)
	            	{
	            		g2.setClip(original_clip);
	            	}
	            }           
	                        
	            if (vis.drawBetaCatenin)
	            {                   
	                int clipx[] = new int[3];
	                int clipy[] = new int[3];
	                
	                for (int node=0; node<3; node++)
	                {
	                	clipx[node] = vertex[node].x;
	                	clipy[node] = vertex[node].y;
	                }
	                Polygon clip = new Polygon(clipx,clipy,3);
	                boolean clip_me = false;
	                
	                // Is circumcentre in the triangle? If not, then we'll clip 
	                // the next bit of drawing to fit inside the triangle 
	                if (!clip.contains(new Point(plotcircumcentre.x, plotcircumcentre.y)))
	                {
	                	clip_me = true;
	                	g2.setClip(clip);
	                } 
	                
	                // Plot membrane-bound beta catenin levels
	                for (int node=0; node<3; node++)
	                {    
	                	 SetCellBetaCateninColour(vis.beta_catenin_values[vis.timeStep][index[node]][0], index[node]);
	                	 int xs[] = new int[4];
	                     int ys[]= new int[4];
	                     xs[0] = plotcircumcentre.x;
	                     ys[0] = plotcircumcentre.y;
	                     xs[1] = midpoint[(node+1)%3].x;
	                     ys[1] = midpoint[(node+1)%3].y;
	                     xs[2] = vertex[node].x;
	                     ys[2] = vertex[node].y;
	                     xs[3] = midpoint[(node+2)%3].x;
	                     ys[3] = midpoint[(node+2)%3].y;
	                     g2.fillPolygon(xs,ys,4);
	                }
	                
	                // Plot cytoplasmic beta catenin levels
	                for (int node=0; node<3; node++)
	                {    
	                	 SetCellBetaCateninColour(vis.beta_catenin_values[vis.timeStep][index[node]][1], index[node]);
	                	 r1 = vis.positions[vis.timeStep][index[0]];
	                	 double cyto_scaler = 0.8;
	                	 double mid_cyto_scaler = (cyto_scaler)/2.0;
	                	 double circumcentre_for_vertex_x = (1-cyto_scaler)*vis.positions[vis.timeStep][index[node]].x+cyto_scaler*circumcentre.x;
	                	 double circumcentre_for_vertex_y = (1-cyto_scaler)*vis.positions[vis.timeStep][index[node]].y+cyto_scaler*circumcentre.y;
	                	 PlotPoint scaled_circumcentre_for_vertex = scale(circumcentre_for_vertex_x, circumcentre_for_vertex_y);
	                	 
	                	 double mid1_for_vertex_x = (1-mid_cyto_scaler)*vis.positions[vis.timeStep][index[node]].x+mid_cyto_scaler*vis.positions[vis.timeStep][index[(node+2)%3]].x;
	                	 double mid1_for_vertex_y = (1-mid_cyto_scaler)*vis.positions[vis.timeStep][index[node]].y+mid_cyto_scaler*vis.positions[vis.timeStep][index[(node+2)%3]].y;
	                	 PlotPoint scaled_mid1_for_vertex = scale(mid1_for_vertex_x, mid1_for_vertex_y);
	                	 
	                	 double mid2_for_vertex_x = (1-mid_cyto_scaler)*vis.positions[vis.timeStep][index[node]].x+mid_cyto_scaler*vis.positions[vis.timeStep][index[(node+1)%3]].x;
	                	 double mid2_for_vertex_y = (1-mid_cyto_scaler)*vis.positions[vis.timeStep][index[node]].y+mid_cyto_scaler*vis.positions[vis.timeStep][index[(node+1)%3]].y;
	                	 PlotPoint scaled_mid2_for_vertex = scale(mid2_for_vertex_x, mid2_for_vertex_y);
	                	 
	                     int xs[] = new int[4];
	                     int ys[] = new int[4];
	                     xs[0] = scaled_circumcentre_for_vertex.x;
	                     ys[0] = scaled_circumcentre_for_vertex.y;
	                     xs[1] = scaled_mid1_for_vertex.x;
	                     ys[1] = scaled_mid1_for_vertex.y;
	                     xs[2] = vertex[node].x;
	                     ys[2] = vertex[node].y;
	                     xs[3] = scaled_mid2_for_vertex.x;
	                     ys[3] = scaled_mid2_for_vertex.y;
	                     g2.fillPolygon(xs,ys,4);
	                }
	                
	                g2.setColor(Color.black);
	                
	                // Plot membrane-bound beta catenin levels
	                if ( (vis.cell_type[vis.timeStep][index[0]] != INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[1]] != INVISIBLE_COLOUR) )
	                {
	                	g2.drawLine(midpoint[2].x, midpoint[2].y, plotcircumcentre.x, plotcircumcentre.y);
	                }
	                if ( (vis.cell_type[vis.timeStep][index[1]] != INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[2]] != INVISIBLE_COLOUR) )
	                {
	                	g2.drawLine(midpoint[0].x, midpoint[0].y, plotcircumcentre.x, plotcircumcentre.y);
	                }
	                if ( (vis.cell_type[vis.timeStep][index[2]] != INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[0]] != INVISIBLE_COLOUR) )
	                {
	                	g2.drawLine(midpoint[1].x, midpoint[1].y, plotcircumcentre.x, plotcircumcentre.y);
	                }
	                if (clip_me)
	                {
	                	g2.setClip(original_clip);
	                }
	            }       
	            
	            if (vis.drawSprings)
	            {
	                // Plot lines
	                if ( (vis.cell_type[vis.timeStep][index[0]] != INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[1]] != INVISIBLE_COLOUR) )
	                {
	                    g2.drawLine(vertex[0].x, vertex[0].y, vertex[1].x, vertex[1].y);
	                }
	                if ( (vis.cell_type[vis.timeStep][index[1]] != INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[2]] != INVISIBLE_COLOUR) )
	                {
	                    g2.drawLine(vertex[1].x, vertex[1].y, vertex[2].x, vertex[2].y);
	                }
	                if ( (vis.cell_type[vis.timeStep][index[2]] != INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[0]] != INVISIBLE_COLOUR) )
	                {
	                    g2.drawLine(vertex[2].x, vertex[2].y, vertex[0].x, vertex[0].y);
	                }
	                if (vis.drawGhosts)
	                {
	                    g2.setColor(spring_silver);
	                    if ( (vis.cell_type[vis.timeStep][index[0]] == INVISIBLE_COLOUR) || (vis.cell_type[vis.timeStep][index[1]] == INVISIBLE_COLOUR) )
	                    {
	                        g2.drawLine(vertex[0].x, vertex[0].y, vertex[1].x, vertex[1].y);
	                    }
	                    if ( (vis.cell_type[vis.timeStep][index[1]] == INVISIBLE_COLOUR) || (vis.cell_type[vis.timeStep][index[2]] == INVISIBLE_COLOUR) )
	                    {
	                        g2.drawLine(vertex[1].x, vertex[1].y, vertex[2].x, vertex[2].y);
	                    }
	                    if ( (vis.cell_type[vis.timeStep][index[2]] == INVISIBLE_COLOUR) || (vis.cell_type[vis.timeStep][index[0]] == INVISIBLE_COLOUR) )
	                    {
	                        g2.drawLine(vertex[2].x, vertex[2].y, vertex[0].x, vertex[0].y);
	                    }
	                    g2.setColor(Color.black);
	                }
	            }
	        }
        }
	        
        if (vis.elementFilePresent && vis.isSparseMesh )
        {        	
	        //  1D Elements in 2D space 
	        for (int i=0; i<vis.numElements[vis.timeStep]; i++)
	        {       
	            // What nodes are we joining up?
	        	int index[] = new int[2];
	            index[0] = vis.element_nodes[vis.timeStep][2*i];
	            index[1] = vis.element_nodes[vis.timeStep][2*i+1];
	            
	            
	            RealPoint r1 = vis.positions[vis.timeStep][index[0]];
	            RealPoint r2 = vis.positions[vis.timeStep][index[1]];
	            
	            // Where are they? Convert to integer pixels
	            PlotPoint vertex[] = new PlotPoint[3];
	            vertex[0] = scale(r1);
	            vertex[1] = scale(r2);
	                        
	            g2.setColor(Color.black);
	            
	            if (vis.drawSprings)
	            {
	                // Plot lines
	                if ( (vis.cell_type[vis.timeStep][index[0]] != INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[1]] != INVISIBLE_COLOUR) )
	                {
	                    g2.drawLine(vertex[0].x, vertex[0].y, vertex[1].x, vertex[1].y);
	                }
	                if (vis.drawGhosts)
	                {
	                    g2.setColor(spring_silver);
	                    if ( (vis.cell_type[vis.timeStep][index[0]] == INVISIBLE_COLOUR) || (vis.cell_type[vis.timeStep][index[1]] == INVISIBLE_COLOUR) )
	                    {
	                        g2.drawLine(vertex[0].x, vertex[0].y, vertex[1].x, vertex[1].y);
	                    }
	                    g2.setColor(Color.black);
	                }
	            }
	        }
        }

        // Draw nodes second so that dots are on top of lines
        double fibre_length = 1.2*node_radius;
        double orientation_length = 5*node_radius;
        for (int i=0; i<vis.numCells[vis.timeStep]; i++) 
        {
        	PlotPoint p = scale(vis.positions[vis.timeStep][i]);

        	if (vis.drawBetaCatenin)
            { 
        		SetCellBetaCateninColour(vis.beta_catenin_values[vis.timeStep][i][2], i);
        	}
        	else
        	{
        		SetNodeColour(i);
        	}
        	if (!vis.drawNutrient)
        	{
        		// \todo: Larger simulations would be clearer with smaller nodes
        		g2.fillOval(p.x - node_radius, p.y - node_radius, 2 * node_radius, 2 * node_radius);
        	}
        	if (vis.drawOrientations)
            {
        		g2.setColor(Color.magenta);
        		
        		RealPoint orientation = vis.orientations[vis.timeStep][i];
        	    g2.drawLine((int) (p.x-orientation_length*orientation.x), (int) (p.y+orientation_length*orientation.y), (int) (p.x+orientation_length*orientation.x), (int) (p.y-orientation_length*orientation.y) );
            }
        	if (vis.drawFibres)
        	{
        		g2.setColor(Color.magenta);
        		RealPoint fibre = vis.fibres[vis.timeStep][i];
        		g2.drawLine(p.x, p.y, (int) (p.x+fibre_length*fibre.x), (int) (p.y-fibre_length*fibre.y) );
        	}
        }
        g2.setColor(Color.black);
        
        if (vis.drawAxes)
        {
        	drawXAxis(tick_length);
        	drawYAxis(tick_length);
        }
        
        if (vis.drawNutrient)
        {
        	drawNutrientColourBar(); 
        }
        if (vis.drawBetaCatenin)
        {
        	drawBetaCateninColourBar(); 
        }
        if (vis.drawAverageStress || vis.drawDifferenceStress)
        {
        	drawStressColourBar(); 
        }
        imageReady = true;
    }
    
    private void drawNutrientColourBar() 
    {
    	int panelHeight = (int) (0.8 * vis.frame.getHeight());
        int panelWidth = (int) (0.8 * vis.frame.getWidth());

        int num_blocks = 100;
        int blockWidth = panelWidth/25;
        int blockHeight = panelHeight/(2*num_blocks);

    	// Draw nutrient colour bar
        for (int i=1; i<=num_blocks; i++)
        {       
            // Calculate colour 
            double conc = (double)(num_blocks - i)/(double)(num_blocks);
            int g = (int)(255.0 * conc);
            int b = (int)(200.0 - 80.0*conc); 

            g2.setColor(new Color(0,g,b));
            g2.fillRect(panelWidth, i*blockHeight+5, blockWidth, blockHeight);
        }

        // Populate vector of labels for nutrient colour bar
        // If num_labels > 11 then change the decimal format as appropriate
        int num_labels = 11;
    	String[] labels = new String[num_labels];
        NumberFormat formatter = new DecimalFormat("#0.0");

    	for (int i=0; i<num_labels; i++)
    	{
    		double this_label_as_number = ((double) num_labels-1-i)*(1.0 / (double) (num_labels-1));
    		String this_label_as_string = formatter.format(this_label_as_number).toString();
    		labels[i] = this_label_as_string;
    	}

        // Draw labels next to nutrient colour bar
    	int block_multiple = num_blocks / (num_labels-1);
        for (int i=0; i<num_labels; i++)
        {
            g2.setColor(Color.black);
            g2.drawString(labels[i], panelWidth + blockWidth + 5, (1 + block_multiple*i)*blockHeight + 5); 
        }     
    }
    
    private void drawBetaCateninColourBar()  
    {
    	int panelHeight = (int) (0.8 * vis.frame.getHeight());
        int panelWidth = (int) (0.8 * vis.frame.getWidth());
        int num_blocks = 10;        
        int blockHeight = panelHeight/(2*num_blocks);
        int blockWidth = blockHeight;
        
        String[] labels = {"20", "18", "16", "14", "12", "10", "8", "6", "4", "2", "0"};        
               
        for (int i=1; i<=num_blocks; i++)
        {       
            // Calculate colour 
        	int g = Math.min(100+16*(num_blocks-i+1),255);
            Color colour = new Color(100,g,100);
            g2.setColor(colour);
            
            g2.fillRect(panelWidth, i*blockHeight, blockWidth, blockHeight);
            g2.setColor(Color.black);
            g2.drawString(labels[i], panelWidth + blockWidth + 5, (i+1)*blockHeight);            
        }   
        g2.drawString(labels[0], panelWidth + blockWidth + 5, blockHeight);  
    }
    
    private void drawStressColourBar()  
    {
    	int panelHeight = (int) (0.8 * vis.frame.getHeight());
        int panelWidth = (int) (0.8 * vis.frame.getWidth());
        int num_blocks = 10;        
        int blockHeight = panelHeight/(2*num_blocks);
        int blockWidth = blockHeight;
        
        double stress_max = 40;
    	double stress_min = 0;
    	int r = 0;
    	int g = 0;
    	int b = 0;
    	double interval = (stress_max-stress_min)/4.0;
        
        for (int i=1; i<=num_blocks; i++)
        {       
            // Calculate colour 
        	double stress =  stress_min + (stress_max-stress_min)*i/(num_blocks-1);        	
            if (stress < interval)
        	{
        		r = 0;
        		g = Math.min( Math.max((int)(255.0*stress/interval),0) , 255 );  
        		b = 255;
        	}
        	else if (stress < 2.0*interval)
        	{
        		r = 0;
        		g = 255;
        		b = Math.min( Math.max(255 - (int) (255.0*(stress-interval)/interval),0) , 255 );
        	}
        	else if (stress < 3.0*interval)
        	{
        		r = Math.min( Math.max((int)(255.0*(stress-2.0*interval)/interval),0) , 255 );
        		g = 255;
        		b = 0;
        	}
        	else 
        	{
        		r = 255;
        		g = Math.min( Math.max(255 - (int) (255.0*(stress- 3.0*interval)/interval),0) , 255 );
        		b = 0;
        	}
            
            Color colour = new Color(r,g,b);
            g2.setColor(colour);
            
            g2.fillRect(panelWidth, i*blockHeight, blockWidth, blockHeight);
            g2.setColor(Color.black);
            
            double box_interval = (stress_max-stress_min)/(double)(num_blocks);
            double lower_bound = stress_min + i*box_interval;
            double upper_bound = lower_bound + box_interval;
            String colour_label = new String(lower_bound+ " - " + upper_bound); 
            g2.drawString(colour_label, panelWidth + blockWidth + 5, (i+1)*blockHeight);            
        }           
    }
    
    private void drawXAxis(int tick_length) 
    {
        int min_x = (int) vis.min_x;
        if (vis.min_x<0)
        {
        	min_x -= 1;
        }
        int max_x = (int) vis.max_x;
        if (vis.max_x>0)
        {
        	max_x += 1;
        }
        
        // work out the number of ticks to use - if it's too big (>15, say)
        // this would make the axis look crowded, so halve num_ticks  
        if ((max_x-min_x)%2 != 0)
        {
        	max_x++;
        }
        int num_ticks = max_x - min_x;        
        int tick_spacing = 1;
        if (num_ticks > 11) 
        {        	
        	num_ticks = num_ticks/2;
        	tick_spacing = 2;
        }
        
        PlotPoint start = scale(min_x, 0);
        PlotPoint end = scale(max_x, 0);
        g2.drawLine(start.x, start.y, end.x, end.y);
          
                
        for (int i = 0; i <= num_ticks; i++) 
        {
            double x = (double) (min_x + tick_spacing*i);
            DecimalFormat df = new DecimalFormat("0.0");
            String x_1dp = df.format(x);
              
            // Tick lines
            PlotPoint posn =  scale(x,0);
            g2.drawLine(posn.x, posn.y, posn.x, posn.y+tick_length);
            g2.drawString(x_1dp, posn.x, posn.y + 2*tick_length);
        }
    }
        
    private void drawYAxis(int tick_length) 
    {
        int min_y = (int) vis.min_y;
        if (vis.min_y<0)
        {
        	min_y -= 1;
        }
        int max_y = (int) vis.max_y;
        if (vis.max_y>0)
        {
        	max_y += 1;
        }        
        
        // work out the number of ticks to use - if it's too big (>15, say)
        // this would make the axis look crowded, so halve num_ticks  
        if ((max_y-min_y)%2 != 0)
        {
        	max_y++;
        }
        int num_ticks = max_y - min_y;        
        int tick_spacing = 1;
        if (num_ticks > 11) 
        {        	
        	num_ticks = num_ticks/2;
        	tick_spacing = 2;
        }
                
        PlotPoint start = scale(0, min_y);
        PlotPoint end = scale(0, max_y);        
        g2.drawLine(start.x, start.y, end.x, end.y);
        
        for (int i=0; i<=num_ticks; i++) 
        {
            double y = (double) (min_y + tick_spacing*i);
            DecimalFormat df = new DecimalFormat("0.0");
            String y_1dp = df.format(y);

            // Tick lines
            PlotPoint posn = scale(0,y);
            g2.drawLine(posn.x-tick_length, posn.y, posn.x, posn.y);
            g2.drawString(y_1dp, posn.x - 4*tick_length, posn.y );
        }
    }
    
    PlotPoint scale(double x, double y)
    {
        // Map min_x to eps and max_x to width-eps (to allow a border)
        int eps = 100;
        int ix = (int) ((x - vis.min_x) * (width-2*eps) /(vis.max_x - vis.min_x) +eps);
        int iy = (int) ((y - vis.min_y) * (height-2*eps) /(vis.max_y - vis.min_y) +eps);
        iy = height - iy; // this is because java is silly and has the y axis going down the screen
        return (new PlotPoint(ix, iy));
    }
    
    PlotPoint scale(RealPoint p) 
    {
        return (scale(p.x, p.y));
    }
    
    RealPoint unscale(PlotPoint p)
    {
        int ix = p.x;
        int iy = height - p.y;
        int eps = 20;
        
        double x = (ix-eps)*(vis.max_x - vis.min_x) / (width-2.0*eps) + vis.min_x;
        double y = (iy-eps)*(vis.max_y - vis.min_y) / (height-2.0*eps) + vis.min_y;
        
        return (new RealPoint(x, y));
    }
    
    public void mouseMoved(MouseEvent e) 
    {
        PlotPoint mouse_position = new PlotPoint(e.getX(), e.getY());
        RealPoint real_position = unscale(mouse_position);
        
        int nearest_index = -1;
        for (int i=0; i<vis.numCells[vis.timeStep]; i++) 
        {
            int sq_dist = SquaredDistance(scale(vis.positions[vis.timeStep][i]), mouse_position);
            if (sq_dist < node_radius*node_radius)
            {
                nearest_index = i;
                break;
            }
        }
        if (nearest_index >= 0)
        {
            RealPoint node_position = vis.positions[vis.timeStep][nearest_index];
            vis.nearest_label.setText("Node "+nearest_index+" is at "+ node_position.x + "  " + node_position.y);
        }
        else
        {
            vis.nearest_label.setText("");
        }
    }
    
    public void mouseDragged(MouseEvent e) 
    {
        // Not used
    }
    
    int SquaredDistance(PlotPoint p0, PlotPoint p1)
    {
        int diffx = p0.x-p1.x;
        int diffy = p0.y-p1.y;
        return diffx*diffx + diffy*diffy;
    }
    
    RealPoint DrawCircumcentre(RealPoint p0, RealPoint p1, RealPoint p2)
    {
        /* To find the coordinates (x_c,y_c) of the circumcentre, we first translate
         * coordinates so that one vertex of the element is at the origin and the 
         * other two vertices are at (X1,Y1) and (X2,Y2). We then solve the linear 
         * system
         * 
         *  ( 2*X1 2*Y1 ) (x_c) = (X1^2 + Y1^2)
         *  ( 2*X2 2*Y2 ) (y_c)   (X2^2 + Y2^2)
         * 
         */ 
        
        double X1 = p1.x - p0.x;
        double Y1 = p1.y - p0.y;
        double X2 = p2.x - p0.x;
        double Y2 = p2.y - p0.y;
        double determinant = X1*Y2 - X2*Y1;
        double RHS1 = (X1*X1 + Y1*Y1)/2.0;
        double RHS2 = (X2*X2 + Y2*Y2)/2.0;
        double x_c = (Y2*RHS1 - Y1*RHS2)/determinant;
        double y_c = (-X2*RHS1 + X1*RHS2)/determinant;
        
        // Translate back to original coordinate system
        x_c += p0.x;
        y_c += p0.y;
        
        return (new RealPoint(x_c,y_c));
    }
     
    void SetNodeColour(int index)
    { 
        if (vis.cell_type[vis.timeStep][index]!=INVISIBLE_COLOUR && 
            vis.drawAncestors && 
            vis.ancestor_values[vis.timeStep][index]!=-1)
        {
            Color ancestor_colour = ancestorColourMap(vis.ancestor_values[vis.timeStep][index]);
            int new_r = 0;
            int new_g = 0;
            int new_b = 0;
            if (ancestor_colour.getRed() - 40 > new_r) new_r = ancestor_colour.getRed() - 40;
            if (ancestor_colour.getGreen() - 40 > new_g) new_g = ancestor_colour.getGreen() - 40;
            if (ancestor_colour.getBlue() - 40 > new_b) new_b = ancestor_colour.getBlue() - 40;
            g2.setColor(new Color(new_r, new_g, new_b));
        }
        else
        {
    		switch (vis.cell_type[vis.timeStep][index]) 
    		{
    	    	case STEM_COLOUR: // stem cell
            	    g2.setColor(Color.green); 
            	    break;
           	 	case TRANSIT_COLOUR: // transit cell
            	    g2.setColor(Color.orange); 
            	    break;
	            case DIFFERENTIATED_COLOUR: // differentiated cell
    	            g2.setColor(Color.red); 
    	            break;
    	        case EARLY_CANCER_COLOUR: // early cancer
    	            g2.setColor(Color.gray); 
    	            break;
    	        case LATE_CANCER_COLOUR:  // late cancer
    	            g2.setColor(Color.black);
    	            break;
    	        case LABELLED_COLOUR: // labelled cell
    	            g2.setColor(Color.blue); 
    	            break;
    	        case APOPTOSIS_COLOUR: // apoptotic cell
    	            g2.setColor(Color.black); 
    	            break;
    	        case INVISIBLE_COLOUR: // sloughed cell
    	            if (!vis.drawGhosts)
    	            {
    	                g2.setColor(background_white);
    	            }
    	            else
    	            {
    	                g2.setColor(Color.lightGray);
    	            } 
    	            break;
				default: 
                	g2.setColor(Color.white);
                	break;
			}  
        }  	
    }
 
     
    /* From C++ 
    enum cell_colours
    {
        STEM_COLOUR,//0
        TRANSIT_COLOUR,//1
        DIFFERENTIATED_COLOUR,//2
        EARLY_CANCER_COLOUR,//3
        LATE_CANCER_COLOUR,//4
        LABELLED_COLOUR,//5
        APOPTOSIS_COLOUR,//6
        INVISIBLE_COLOUR, // visualizer treats '7' as invisible
    };
    */

    public static final int STEM_COLOUR = 0;
    public static final int TRANSIT_COLOUR = 1;
    public static final int DIFFERENTIATED_COLOUR = 2;
    public static final int EARLY_CANCER_COLOUR = 3;
    public static final int LATE_CANCER_COLOUR = 4;
    public static final int LABELLED_COLOUR = 5;
    public static final int APOPTOSIS_COLOUR = 6;
    public static final int INVISIBLE_COLOUR = 7;
    
    void SetCellColour(int index)
    {
    	if (vis.drawNutrient)
    	{
    		double conc = vis.nutrient_values[vis.timeStep][index];
    		
    		switch (vis.cell_type[vis.timeStep][index]) 
        	{
        		case APOPTOSIS_COLOUR: // apoptotic cell
        			g2.setColor(apoptotic_grey); 
        			break;
        		case INVISIBLE_COLOUR: // sloughed cell
        			g2.setColor(background_white); 
        			break;
        		case LABELLED_COLOUR: // labelled cell
        			g2.setColor(Color.red);
        			break;	
        		default: // any other cell type
                	int r = 0;
                	int g = Math.min( Math.max((int)(255*conc),0) , 255 );            	
                    int b = Math.min( Math.max( (int)(200 - 80*conc),0) , 255);                  
                    Color colour = new Color(r,g,b);
                    g2.setColor(colour);
        		    break;
        	}
    		
    	}
    	else if (vis.drawAverageStress || vis.drawDifferenceStress)
    	{
    		double stress;
        	if (vis.drawAverageStress)
        	{
        		stress = vis.stress_values[vis.timeStep][index][0];
        	}
        	else
        	{
        		stress = vis.stress_values[vis.timeStep][index][1];
            }
        	
    		double stress_max = 40;
        	double stress_min = 0;
        	int r = 0;
        	int g = 0;
        	int b = 0;
        	double interval = (stress_max-stress_min)/4.0;
        	
        	// We do not show up ghost nodes in a stress plot.
            if (vis.cell_type[vis.timeStep][index] == INVISIBLE_COLOUR)
            {
            	// Sloughed cell
                g2.setColor(background_white);
            }
            else
            {
            	if (stress < interval)
            	{
            		r = 0;
            		g = Math.min( Math.max((int)(255.0*stress/interval),0) , 255 ); 
            		b = 255;
            	}
            	else if (stress < 2.0*interval)
            	{
            		r = 0;
            		g = 255;
            		b = Math.min( Math.max(255 - (int) (255.0*(stress - interval)/interval),0) , 255 ); 
            	}
            	else if (stress < 3.0*interval)
            	{
            		r = Math.min( Math.max((int)(255.0*(stress - 2.0*interval)/interval),0) , 255 ); 
            		g = 255;
            		b = 0;
            	}
            	else 
            	{
            		r = 255;
            		g = Math.min( Math.max(255 - (int) (255.0*(stress - 3.0*interval)/interval),0) , 255 ); 
            		b = 0;
            	}
                Color colour = new Color(r,g,b);
                g2.setColor(colour);
            }   	
    	}
    	else if (vis.drawAncestors && (vis.ancestor_values[vis.timeStep][index]!=-1))
      	{	// If we are drawing ancestors and this cell's value has been set in simulation.
    		if (vis.cell_type[vis.timeStep][index] == INVISIBLE_COLOUR)
    		{
    			g2.setColor(background_white);      			
    		}
    		else
    		{
                        Color ancestor_colour = ancestorColourMap(vis.ancestor_values[vis.timeStep][index]);
    			g2.setColor(ancestor_colour);
    		}
      	}
    	else
        {
    		switch (vis.cell_type[vis.timeStep][index]) 
        	{
        		case STEM_COLOUR: // stem cell
        			g2.setColor(Color.cyan); 
        			break;
        		case TRANSIT_COLOUR: // transit cell
        			g2.setColor(Color.yellow); 
        			break;
        		case DIFFERENTIATED_COLOUR: // differentiated cell
        			g2.setColor(Color.pink); 
        			break;
        		case EARLY_CANCER_COLOUR: // early cancer
        			g2.setColor(Color.lightGray); 
        			break;
        		case LATE_CANCER_COLOUR:  // late cancer
        			g2.setColor(Color.gray);
        			break;
        		case LABELLED_COLOUR: // labelled cell
        			g2.setColor(purple); 
        			break;
        		case APOPTOSIS_COLOUR: // apoptotic cell
        			g2.setColor(apoptotic_grey); 
        			break;
        		case INVISIBLE_COLOUR: // sloughed cell
        			g2.setColor(background_white);                    
        			break;
        		default: 
        			g2.setColor(background_white);                    
        			break;
        	}
    	}
    }

    public Color ancestorColourMap(int ancestor)
    {
        //Map the colour uniquely into [0, 255]
        int r=hash32shiftmult(ancestor, 256);
        int g=hash32shiftmult(ancestor+1, 256);
        int b=hash32shiftmult(ancestor*2, 256);
        return new Color(r,g,b);
    }

    public int hash32shiftmult(int key, int range)
    {
      //Mostly copied from http://www.acme.com/resources/classes/Acme/IntHashtable.java
      int c2=0x27d4eb2d; // a prime or an odd constant
      key = (key ^ 61) ^ (key >>> 16);
      key = key + (key << 3);
      key = key ^ (key >>> 4);
      key = key * c2;
      key = key ^ (key >>> 15);
      //We added the last two lines
      key = key & 0x7FFFFFFF; //Make positive unsigned
      return (key%range);//In 0<=key<range
    }

        
    void SetCellBetaCateninColour(double conc, int index)
    {
    	if (vis.cell_type[vis.timeStep][index] == INVISIBLE_COLOUR)
        {
            // Sloughed cell
            g2.setColor(background_white);
        }
    	else
    	{
    		int r = 100;
    		int b = 100;
    		int g = Math.min( Math.max((int)(100 + 8*conc),0) , 255 );  
    		Color colour = new Color(r,g,b);
    		g2.setColor(colour);
    	}
    }
    
 }
