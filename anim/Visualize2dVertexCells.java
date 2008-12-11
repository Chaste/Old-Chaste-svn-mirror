/*

Copyright (C) University of Oxford, 2008

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
import java.text.DecimalFormat;
import java.util.StringTokenizer;
import java.lang.Math;
import javax.imageio.ImageIO;

import javax.swing.JPanel;
import javax.swing.JLabel;

public class Visualize2dVertexCells implements ActionListener, AdjustmentListener, ItemListener, Runnable
{
    private Thread updateThread;

    static CustomCanvas2D canvas;

    Button run;
    
    public Frame frame = new Frame();
    
    public static boolean parsed_all_files = false;
    public static boolean drawAncestors = false;
    public static boolean drawAxes = true;
    public static boolean drawCells = true;
    public static boolean writeFiles = false;
    public static boolean drawFibres = false;
    public static boolean drawCylinder = false;
    public static boolean drawCylinderOverride = true;
    public static boolean ancestorsFilePresent = false;
    public static boolean elementFilePresent = true;
    // by default the last timestep isn't read or visualised; this
    // allows the visualiser to be run as a simulation is being run.
    // To visualise the last timestep, use "showlaststep" as an argument
    public static boolean showLastStep = false; 
    
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
    public static double[] times;
    
    public static RealPoint[][] positions;
    
    public static Scrollbar delay_slider = new Scrollbar(Scrollbar.HORIZONTAL, delay, 1, 1, 100);
    public static Scrollbar time_slider = new Scrollbar(Scrollbar.HORIZONTAL, timeStep, 1, 0, 2);
    
    public static Checkbox output = new Checkbox("Output");
    public static Checkbox cells = new Checkbox("Cells");
    public static Checkbox ancestors = new Checkbox("Clonal Populations");
    public static Checkbox axes = new Checkbox("Axes");
    
    public static JLabel nearest_label = new JLabel();
    
    public static File node_file;
    public static File element_file;
    public static File ancestors_file;
    public static File setup_file;
    
    public static Button refresh;
    
    public Visualize2dVertexCells() 
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
        else if (cb == cells)
        {
            System.out.println("Drawing cells = "+drawCells);
        }
        else if (cb == axes)
        {
            drawAxes = state;
            System.out.println("Drawing axes = "+drawAxes); 
        }
        else if (cb == ancestors)
        {
            drawAncestors = state;
            System.out.println("Drawing clonal populations = " + drawAncestors); 
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
                    Thread fred = updateThread;
                    updateThread = null;
                    fred.interrupt();
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
        cells.addItemListener(this);
        axes.addItemListener(this);
        ancestors.addItemListener(this);
        
        checkPanel.add(output);
        checkPanel.add(cells);
        checkPanel.add(axes);
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
        cells.setState(true);
        output.setState(false);
        axes.setState(true);
        ancestors.setState(false);
        
        for (int i=1; i<args.length; i++)
        {
            if (args[i].equals("output"))
            {
                writeFiles = true;
                output.setState(true);                  
            } 
            else if (args[i].equals("nocells"))
            {
                drawCells = false;
                cells.setState(false);
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
            else if (args[i].equals("ancestors"))
            {
                drawAncestors = true;
                ancestors.setState(true);
            }
            else
            {
                System.out.println("Input option not recognised");
            }
        }
        node_file = new File(args[0]+"/results.viznodes");
        element_file = new File(args[0]+"/results.vizelements");
        ancestors_file = new File(args[0]+"/results.vizAncestors");
                
        if (!node_file.isFile())
        {
            System.out.println("The file "+args[0]+"/results.viznodes doesn't exist");
            return;
        }
        if (!element_file.isFile())
        {
        	// If the results.vizelements does not exist, then assume the results 
        	// were generated using a SimpleTissue
        	elementFilePresent = false;
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
                        drawCylinder = true && drawCylinderOverride;    // this is made true only if mesh width exists
                    }
                    if (parameter.equals("Complete")) 
                    {
                	    showLastStep = true;
                    }
                    line_setup = in_setup_file.readLine();
                }
            }
            catch (Exception e) 
        	{
            	System.out.println("Error occured. Exception message: "+e.getMessage());
        	}
        }             
        
        final Visualize2dVertexCells vis = new Visualize2dVertexCells();

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
            
            // by default we don't read or print the final line, so
            // the visualiser can be run as a simulation is being run,
            // and the visualiser will work on incomplete data.
            if(num_lines>1 && !showLastStep)
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
            if (ancestorsFilePresent)
            {
            	ancestor_values = new int[num_lines][];
            }
            
            String line_ancestors = "";
            BufferedReader in_ancestors_file = null;
            if (drawAncestors)
            {
               	ancestor_values = new int[num_lines][]; 
               	in_ancestors_file = new BufferedReader(new FileReader(ancestors_file));
               	line_ancestors = in_ancestors_file.readLine();
            }

            BufferedReader in_node_file = new BufferedReader(new FileReader(node_file));
            String line_node = in_node_file.readLine(); // from console input example
            
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
            int row = 0;
            while (line_node != null && row<num_lines) 
            {
            	// Create a StringTokenizer with a colon sign as a delimiter
                StringTokenizer st_node = new StringTokenizer(line_node);
                StringTokenizer st_element = null;
                if (elementFilePresent)
                {
                	st_element = new StringTokenizer(line_element);
                }
                StringTokenizer st_ancestors = null;
                
                Double time = Double.valueOf(st_node.nextToken());
                
                
                if (drawAncestors)
                {
                    st_ancestors = new StringTokenizer(line_ancestors);
                    Double ancestors_time = Double.valueOf(st_ancestors.nextToken());
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
// old stuff
                // Count the number of entries in the node file and check correct 
                int entries = st_node.countTokens();
                                
                if (entries%3 != 0)
                {
                	System.out.println("Warning: Results from time "+time.doubleValue()+" will not be plotted as the corresponding line of the node file is not of the required form: time,x,y,type,x,y,type...");
                	break;
                }
                numCells[row] = entries/3; 

//                if (elementFilePresent)
//                {
//                    // Count the number of entries in the element file and check correct 
//                    entries = st_element.countTokens();
//                    if (entries%3 != 0)
//                    {
//                    	System.out.println("Warning: Results from time "+time.doubleValue()+" will not be plotted as the corresponding line of the element file is not of the required form: time,n1,n2,n3,n1,n2,n3..");
//                    	break;
//                    }
//                    numElements[row] = entries/3;
//                    element_nodes[row] = new int[memory_factor*3*numElements[row]];
//                } 
                
                // new stuff
                if (elementFilePresent)
                {
                    // Count the number of entries in the element file and check correct 
                	int total_entries = st_element.countTokens();
                    int total_elems = 0;
                    int entry_posn = 0;
                    System.out.println("About to loop over elems");
                	
                    while (entry_posn<total_entries)
                    {
                     	entry_posn ++;
                     	total_elems++;
                    	int num_elem_vertices = Integer.parseInt(st_element.nextToken());
                    	for (int i=0; i<num_elem_vertices; i++)
                    	{
                    		st_element.nextToken();
                         	entry_posn ++;
                    	}
                    }
                    
                    new int[total_elems];
                    
                    System.out.println("num elements = " + total_elems);
                    numElements[row] = total_elems;
                    element_nodes[row] = new int[memory_factor*(total_entries-total_elems)];
                }                

                positions[row] = new RealPoint[memory_factor*numCells[row]];
                if (ancestorsFilePresent)
                {
                	ancestor_values[row] = new int[memory_factor*numCells[row]];
                }
                
                cell_type[row] = new int[memory_factor*numCells[row]];
                        
                for (int i=0; i<numCells[row]; i++) 
                {
                    double d1 = Double.valueOf(st_node.nextToken()).doubleValue();
                    double d2 = Double.valueOf(st_node.nextToken()).doubleValue();

                    cell_type[row][i] = Integer.parseInt(st_node.nextToken());
                    if ( cell_type[row][i] < 0 )
                    {
                        System.out.println("Error: Cell type must be a non-negative integer");
                        System.exit(0);
                    }
                    positions[row][i] = new RealPoint(d1,d2);
                    
                    if (drawAncestors)
                    {	// If this is a real cell then read in ancestor from row
                    	if (cell_type[row][i]!=canvas.INVISIBLE_COLOUR)	// if this is not a ghost cell
                    	{
                    		int ancestor_value = Integer.parseInt(st_ancestors.nextToken()); // index
                        	ancestor_values[row][i] = ancestor_value;
                        }
                    }	
                    
                    
                    
                }                
                
                if (elementFilePresent)
                {
                	for (int i = 0; i < 3*numElements[row]; i++) 
                    {
                        int node = Integer.parseInt(st_element.nextToken());
                        element_nodes[row][i] = node;
                    }                	
                }
                
                // Read next line of the file
                line_node = in_node_file.readLine();
                
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
            System.out.println("Drawing cells = "+drawCells);
            System.out.println("Drawing axes = "+ drawAxes);
            System.out.println("Drawing clonal populations = "+ drawAncestors);
            
            
            
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
        for (int row=0 ; row<numSteps; row++)
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
        
        //If node C is far away 
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

    Visualize2dVertexCells vis;

    boolean imageReady = false;   
    boolean imageDrawing = false;
    
    int width;
    int height;
    int node_radius = 5;
    
    BufferedImage buffered_image = null;
    Graphics g2 = null;
    
    Color garysSexySilver = new Color(238,238,238);
    Color garysSpringsSilver = new Color(200,200,200);
    Color ozzysDirtyGrey = new Color(80,80,80);
    Color purple = new Color(121,126,234);
    
    public CustomCanvas2D(Visualize2dVertexCells v) 
    {
        vis = v;
        setBackground(garysSexySilver);
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
        
        g2.setColor(garysSexySilver);
        g2.fillRect(0, 0, width, height);
        g2.setColor(Color.black);
        g2.drawString("Time = " + vis.times[vis.timeStep], 10, 10);

        g2.setColor(Color.black);
        Shape original_clip = g2.getClip();
        
        if (vis.elementFilePresent)
        {        	
	        // Draw elements first
	        for (int i=0 ; i < vis.numElements[vis.timeStep]; i++)
	        {       
	            // What nodes are we joining up?
	        	int index[] = new int[3];
	            index[0] = vis.element_nodes[vis.timeStep][3*i];
	            index[1] = vis.element_nodes[vis.timeStep][3*i+1];
	            index[2] = vis.element_nodes[vis.timeStep][3*i+2];
	            
	            RealPoint r1 = vis.positions[vis.timeStep][index[0]];
	            RealPoint r2 = vis.positions[vis.timeStep][index[1]];
	            RealPoint r3 = vis.positions[vis.timeStep][index[2]];
	            
	            RealPoint circumcentre=DrawCircumcentre(r1,r2,r3);
	            PlotPoint plotcircumcentre = scale(circumcentre);
	            
	            // Where are they? Convert to integer pixels
	            PlotPoint vertex[] = new PlotPoint[3];
	            vertex[0] = scale(r1);
	            vertex[1] = scale(r2);
	            vertex[2] = scale(r3);
	
	            PlotPoint midpoint[] = new PlotPoint[3];
	            midpoint[2] = scale(new RealPoint(r1,r2));
	            midpoint[0] = scale(new RealPoint(r2,r3));
	            midpoint[1] = scale(new RealPoint(r3,r1));
	            
	            g2.setColor(Color.black);
	            
//	            if (vis.drawCells)
//	            {
//	                int clipx[] = new int[3];
//	                int clipy[] = new int[3];
//	                for (int node=0; node<3; node++)
//	                {
//	                	clipx[node] = vertex[node].x;
//	                    clipy[node] = vertex[node].y;
//	                }
//	                Polygon clip = new Polygon(clipx,clipy,3);
//	                boolean clip_me = false;
//	                 
//	                // Is circumcentre in the triangle? If not, then we'll clip 
//	                // the next bit of drawing to fit inside the triangle 
//	                if (!clip.contains(new Point(plotcircumcentre.x, plotcircumcentre.y)))
//	                {
//	                	clip_me = true;
//	                    g2.setClip(clip);
//	                }
//	                for (int node=0; node<3; node++)
//	                {                	 
//	                    SetCellColour(index[node]);
//	                    int xs[] = new int[4];
//	                    int ys[] = new int[4];
//	                    xs[0] = plotcircumcentre.x;
//	                    ys[0] = plotcircumcentre.y;
//	                    xs[1] = midpoint[(node+1)%3].x;
//	                    ys[1] = midpoint[(node+1)%3].y;
//	                    xs[2] = vertex[node].x;
//	                    ys[2] = vertex[node].y;
//	                    xs[3] = midpoint[(node+2)%3].x;
//	                    ys[3] = midpoint[(node+2)%3].y;
//	                    g2.fillPolygon(xs,ys,4);
//	                }
//	                
//	                g2.setColor(Color.black);
//	                
//	                // Plot cell boundary lines
//	                if ( (vis.cell_type[vis.timeStep][index[0]]!= INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[1]]!= INVISIBLE_COLOUR) )
//	                {
//	                    g2.drawLine(midpoint[2].x, midpoint[2].y, plotcircumcentre.x, plotcircumcentre.y);
//	                }
//	                if ( (vis.cell_type[vis.timeStep][index[1]]!= INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[2]]!= INVISIBLE_COLOUR) )
//	                {
//	                    g2.drawLine(midpoint[0].x, midpoint[0].y, plotcircumcentre.x, plotcircumcentre.y);
//	                }
//	                if ( (vis.cell_type[vis.timeStep][index[2]]!= INVISIBLE_COLOUR) && (vis.cell_type[vis.timeStep][index[0]]!= INVISIBLE_COLOUR) )
//	                {
//	                    g2.drawLine(midpoint[1].x, midpoint[1].y, plotcircumcentre.x, plotcircumcentre.y);
//	                }
//	                if (clip_me)
//	                {
//	                	g2.setClip(original_clip);
//	                }                
//	            }	        
            
	            
	            if (vis.drawCells)
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
	                
	            }
	        }
        }

        // Draw nodes second so that dots are on top of lines
        for (int i = 0; i < vis.numCells[vis.timeStep]; i++ ) 
        {
        	PlotPoint p = scale(vis.positions[vis.timeStep][i]);

        	SetNodeColour(i);
        }
        g2.setColor(Color.black);
        
        if (vis.drawAxes)
        {
        	drawXAxis(tick_length);
        	drawYAxis(tick_length);
        }
        
        imageReady = true;
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
        
        for (int i = 0; i <= num_ticks; i++) 
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
        return (new PlotPoint(ix,iy));
    }
    
    PlotPoint scale(RealPoint p) 
    {
        return (scale(p.x,p.y));    
    }
    
    RealPoint unscale(PlotPoint p)
    {
        int ix = p.x;
        int iy = height - p.y;
        int eps = 20;
        
        double x = (ix-eps)*(vis.max_x - vis.min_x) / (width-2.0*eps) + vis.min_x;
        double y = (iy-eps)*(vis.max_y - vis.min_y) / (height-2.0*eps) + vis.min_y;
        
        return (new RealPoint(x,y));
    }
    
    public void mouseMoved(MouseEvent e) 
    {
        PlotPoint mouse_position = new PlotPoint(e.getX(), e.getY());
        RealPoint real_position = unscale(mouse_position);
        
        int nearest_index = -1;
        for (int i = 0; i < vis.numCells[vis.timeStep]; i++) 
        {
            int sq_dist = SquaredDistance(scale(vis.positions[vis.timeStep][i]),mouse_position);
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
                g2.setColor(Color.lightGray);
                break;
             default: 
                g2.setColor(Color.white);
                break;
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
    	
    	if (vis.drawAncestors && (vis.ancestor_values[vis.timeStep][index]!=-1))
      	{	// If we are drawing ancestors and this cell's value has been set in simulation.
    		if (vis.cell_type[vis.timeStep][index] == INVISIBLE_COLOUR)
    		{
    			g2.setColor(garysSexySilver);      			
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
        			g2.setColor(ozzysDirtyGrey); 
        			break;
        		case INVISIBLE_COLOUR: // sloughed cell
        			g2.setColor(garysSexySilver);                    
        			break;
        		default: 
        			g2.setColor(garysSexySilver);                    
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

        
   
 }
