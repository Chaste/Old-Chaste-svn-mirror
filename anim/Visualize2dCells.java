import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;
import java.lang.Math;
import javax.imageio.ImageIO;

import javax.swing.Box;
import javax.swing.JPanel;
import javax.swing.JLabel;

public class Visualize2dCells implements ActionListener, AdjustmentListener, ItemListener, Runnable{

    public Frame frame = new Frame();

    public static double[] times;

    public static int[] numCells;
    public static int[] numElements;
    public static int memory_factor = 2;

    public static RealPoint[][] positions;
    public static RealPoint[][] fibres;
    public static double[][][] beta_catenin_values;
    public static int[][] element_nodes;
    public static int[][] cell_type;
    public static int[][] image_cells;
        
    public static double max_x = -1e12;
    public static double max_y = -1e12;
    public static double min_x =  1e12;
    public static double min_y =  1e12;
    public static double crypt_width = 0.0;
    public static double half_width = 0.0;
    
    public static boolean parsed_all_files=false;
    
    public static boolean drawSprings = false;
    public static boolean drawCells = true;
    public static boolean drawCircles = false;
    public static boolean drawNutrient = false;
    public static boolean drawBetaCatenin = false;
    public static boolean writeFiles = false;
    public static boolean drawGhosts = false;
    public static boolean drawFibres = false;
    public static boolean drawCylinder = false;
    public static boolean drawCylinderOverride = true;
    public static boolean setupFilePresent = false;
    
    public static int timeStep = 0;

    public static int delay = 50;
    
    private Thread updateThread;

    static CustomCanvas2D canvas;

    Button run;
    
    public static Scrollbar delay_slider = new Scrollbar(Scrollbar.HORIZONTAL, delay, 1, 1, 100);
    public static Scrollbar time_slider = new Scrollbar(Scrollbar.HORIZONTAL, timeStep, 1, 0, 2);
    public static JPanel nutrient_colour_bar = new JPanel();
    public static JPanel beta_catenin_colour_bar = new JPanel();
    //public static Checkbox output, springs, fibre, cells, ghost_nodes;
    public static Checkbox output = new Checkbox("Output");
    public static Checkbox springs = new Checkbox("Springs");
    public static Checkbox fibre = new Checkbox("Fibres");
    public static Checkbox cells = new Checkbox("Cells");
    public static Checkbox ghost_nodes = new Checkbox("Ghosts");
    public static Checkbox circles = new Checkbox("Cells as circles");
    public static Checkbox nutrient = new Checkbox("Nutrient");
    public static Checkbox beta_catenin = new Checkbox("Beta catenin");
    public static JLabel nearest_label = new JLabel();
    public static int numSteps = 0;

    public static String nutrient_file;
        
    public Visualize2dCells() {
        frame.setSize(700, 700);
        frame.setTitle("Gavaghan's goons' visualization tools (TM)");

        frame.setLayout(new BorderLayout());
        canvas = new CustomCanvas2D(this);
        //canvas.setSize(1000, 400);
        canvas.setPreferredSize(new Dimension(frame.getWidth(),frame.getHeight()));
        canvas.addMouseMotionListener(canvas);
        addButtons(frame);
        addTimeSlider(frame);
        addNutrientColourBar(frame);
        addBetaCateninColourBar(frame);
        JPanel canvasPanel = new JPanel();
        canvasPanel.add(canvas);
        frame.add(canvasPanel, BorderLayout.CENTER);

        frame.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                System.exit(0);
            }
        });
        frame.pack();
        frame.setVisible(true);
    }

    public void actionPerformed(ActionEvent event) {
        String pressed = event.getActionCommand();
        if (pressed == "Quit") {
            frame.dispose();
        }
        if (pressed == "Run") {
            if (timeStep == numSteps - 1) {
                timeStep = 0;
                time_slider.setValue(timeStep);
            }
            if (updateThread == null) {
                run.setLabel("Pause");
                updateThread = new Thread(this);
                updateThread.start();
            }
        }
        if (pressed == "Reset") {
            timeStep = 0;
            time_slider.setValue(timeStep);
            canvas.drawBufferedImage();
            canvas.repaint();
        }
        if (pressed == "Pause") {
            if (updateThread != null) {
                Thread fred = updateThread;
                updateThread = null;
                fred.interrupt();
                run.setLabel("Run");
            }
        }
    }
    public void itemStateChanged(ItemEvent e) {
        Object cb= e.getItemSelectable();
        boolean state=(e.getStateChange() == ItemEvent.SELECTED);
        
        if (cb == output) 
        {
            writeFiles=state;
            System.out.println("Writing output files = "+writeFiles);
        } 
        else if (cb == springs) 
        {
            drawSprings=state;
            System.out.println("Drawing springs = "+drawSprings);
        }
        else if (cb == fibre)
        {
            drawFibres=state;
            System.out.println("Drawing fibres = "+drawFibres);
        }
        else if (cb == cells)
        {
            drawCells=state;
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
            System.out.println("Drawing cells as circles = "+drawCircles);    
        }
        else if (cb == nutrient)
        {
            drawNutrient = state;
            System.out.println("Drawing nutrient = "+drawNutrient); 
        }
        else if (cb == beta_catenin)
        {
            drawBetaCatenin = state;
            System.out.println("Drawing beta catenin = "+drawBetaCatenin); 
        }
        canvas.drawBufferedImage();
        canvas.repaint();
    }

    public void adjustmentValueChanged(AdjustmentEvent e) {
        delay = delay_slider.getValue();
        timeStep = time_slider.getValue();
        //System.out.println("delay = "+delay);
        //System.out.println("timeStep = "+timeStep);
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
        JPanel buttonPanel = new JPanel(new GridLayout(0,3));
        Button quit = new Button("Quit");
        quit.addActionListener(this);

        run = new Button("Run");
        run.addActionListener(this);
        
        Button reset = new Button("Reset");
        reset.addActionListener(this);
        
        buttonPanel.add(quit);
        buttonPanel.add(run);
        buttonPanel.add(reset);
                
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
                        
        JPanel checkPanel = new JPanel(new GridLayout(0,4));
        output.addItemListener(this);
        springs.addItemListener(this);
        fibre.addItemListener(this);
        cells.addItemListener(this);
        ghost_nodes.addItemListener(this);
        circles.addItemListener(this);
        nutrient.addItemListener(this);
        beta_catenin.addItemListener(this);
        
        checkPanel.add(output);
        checkPanel.add(springs);
        checkPanel.add(nutrient);
        checkPanel.add(beta_catenin);
        checkPanel.add(fibre);
        checkPanel.add(cells);
        checkPanel.add(ghost_nodes);
        checkPanel.add(circles);
        
        checkPanel.add(nearest_label);
                
        JPanel southPanel = new JPanel(new GridLayout(2,0));
        
        southPanel.add(scrollPanel_time);
        southPanel.add(checkPanel);
        frame.add(southPanel,BorderLayout.SOUTH);
    }
    
    
    public void addNutrientColourBar(Frame frame) 
    {
        String[] labels = {"0.9-1", "0.8-0.9", "0.7-0.8", "0.6-0.7", "0.5-0.6", "0.4-0.5", "0.3-0.4", "0.2-0.3", "0.1-0.2", "0.0-0.1"};
        
        int panelHeight = (int) (0.8 * frame.getHeight());
        int panelWidth = 120;
        int num_blocks = 10;        
        int blockHeight = panelHeight/num_blocks;
        
        nutrient_colour_bar.setVisible(false);
        nutrient_colour_bar.setPreferredSize(new Dimension(panelWidth,panelHeight));
        nutrient_colour_bar.setLayout(new GridLayout(10,2));
       
        for (int i=num_blocks-1;i>=0;i--)
        {       
            JPanel colour_block = new JPanel();
            colour_block.setPreferredSize(new Dimension(panelWidth/2,panelHeight/num_blocks));
            
            //Calculate colour                    
            colour_block.setBackground(new Color(121-8*i,126-8*i,200-8*i));
            
            Label colour_label = new Label(labels[9-i]);        
        
            nutrient_colour_bar.add(colour_block);
            nutrient_colour_bar.add(colour_label);      
        }
                
        JPanel eastPanel = new JPanel(new GridLayout(1,1));
        
        eastPanel.add(nutrient_colour_bar);
        frame.add(eastPanel,BorderLayout.EAST);
    }
    
    public void addBetaCateninColourBar(Frame frame) 
    {
        String[] labels = {"18-20", "16-18", "14-16", "12-14", "10-12", "8-10", "6-8", "4-6", "2-4", "0-2"};        
        
        int panelHeight = (int) (0.8 * frame.getHeight());
        int panelWidth = 120;        
        int num_blocks = 10;        
        int blockHeight = panelHeight/num_blocks;
        
        beta_catenin_colour_bar.setVisible(false);
        beta_catenin_colour_bar.setPreferredSize(new Dimension(panelWidth,panelHeight));
        beta_catenin_colour_bar.setLayout(new GridLayout(10,2));
       
        for (int i=num_blocks-1;i>=0;i--)
        {       
            JPanel colour_block = new JPanel();
            colour_block.setPreferredSize(new Dimension(panelWidth/2,panelHeight/num_blocks));
            
            //Calculate colour                    
            colour_block.setBackground(new Color(100,100+16*i,100));
            
            Label colour_label = new Label(labels[9-i]);        
        
            beta_catenin_colour_bar.add(colour_block);
            beta_catenin_colour_bar.add(colour_label);      
        }
                
        JPanel westPanel = new JPanel(new GridLayout(1,1));
        
        westPanel.add(beta_catenin_colour_bar);
        frame.add(westPanel,BorderLayout.WEST);
    }
    

    public static void main(String args[]) 
    {
        System.out.println("Copyright Gavaghan's goons");//(Gary Mirams, Sarah Eastburn, Pras Pathmanathan, Alex Fletcher & Joe Pitt-Francis)");
        output.setState(false);
        springs.setState(false);
        fibre.setState(false);
        cells.setState(true);
        ghost_nodes.setState(false);
        circles.setState(false);
        nutrient.setState(false);
        beta_catenin.setState(false);
        
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
            }
            else if (args[i].equals("betacatenin"))
            {
            	drawBetaCatenin = true;
                beta_catenin.setState(true);
            }
            else if (args[i].equals("notcylindrical"))
            {
                drawCylinderOverride = false;
            }
            else
            {
                System.out.println("Input option not recognised");
            }
        }
        File node_file = new File(args[0]+"/vis_results/results.viznodes");
        File element_file = new File(args[0]+"/vis_results/results.vizelements");
        File beta_catenin_file = new File(args[0]+"/vis_results/results.vizbCat");
        
        // save where the nutrient file will be, in case it is needed later
        nutrient_file = args[0] + "/../nutrients/";

                
        if (!node_file.isFile())
        {
            System.out.println("The file "+args[0]+"/vis_results/results.viznodes doesn't exist");
            return;
        }
        if (!element_file.isFile())
        {
            System.out.println("The file "+args[0]+"/vis_results/results.vizelements doesn't exist");
            return;
        }
        if (!beta_catenin_file.isFile())
        {
            System.out.println("The file "+args[0]+"/vis_results/results.vizbCat doesn't exist");
            return;
        }
    
        File fibre_file= new File(args[0]+"/vis_results/results.vizfibres");
        if (!fibre_file.isFile())
        {
            System.out.println("The file "+args[0]+"/vis_results/results.vizfibres doesn't exist");
            fibre.setVisible(false);
            drawFibres=false;
        } 
        else 
        {
            fibre.setState(true);
            drawFibres=true; //Sorry, this is just to get it working
        }
        
        File setup_file = new File(args[0]+"/vis_results/results.vizsetup");
        if (!setup_file.isFile())
        {
            System.out.println("The file "+args[0]+"/vis_results/results.vizsetup doesn't exist");
        }
        else 
        {
            setupFilePresent = true;
        }
        
        Visualize2dCells vis = new Visualize2dCells();
        
        try 
        {
            BufferedReader skim_node_file = new BufferedReader(new FileReader(node_file));

            int num_lines = 0;
            while (skim_node_file.readLine() != null) 
            {
                num_lines++;
            }

            numSteps = num_lines;
            time_slider.setMaximum(numSteps);
            times = new double[num_lines];
            positions = new RealPoint[num_lines][];
            cell_type = new int [num_lines][];
            numCells = new int[num_lines];
            numElements = new int[num_lines];
            element_nodes = new int[num_lines][];
            image_cells = new int[num_lines][];            
            fibres =  new RealPoint[num_lines][];
            beta_catenin_values = new double[num_lines][][]; 
            String line_fibre="";
            BufferedReader in_fibre_file=null;
            if (drawFibres)
            {
                fibres = new RealPoint[num_lines][]; 
                in_fibre_file=new BufferedReader(new FileReader(fibre_file));
                line_fibre=in_fibre_file.readLine();
            }
                        
            if (setupFilePresent)
            {
                BufferedReader in_setup_file = new BufferedReader(new FileReader(setup_file));
                String line_setup = in_setup_file.readLine();   // above.
                // Read setup information.
                while (line_setup != null)
                {
                    StringTokenizer st_setup = new StringTokenizer(line_setup);
                    String parameter = st_setup.nextToken();
                    if (parameter.equals("MeshWidth"))  // .equals?? That took some doing!
                    {
                        crypt_width = Double.valueOf(st_setup.nextToken());
                        half_width = crypt_width/2.0;
                        System.out.println("Mesh Width = " + crypt_width);
                        drawCylinder = true && drawCylinderOverride;    // this is made true only if mesh width exists. 
                    }
                    if (parameter.equals("BetaCatenin"))  // .equals?? That took some doing!
                    {
                        System.out.println("Using Beta Catenin");
                        drawBetaCatenin = true;
                        beta_catenin.setState(true);
                    }
                    line_setup = in_setup_file.readLine();
                }
            }

            String line_beta_catenin="";
            BufferedReader in_beta_catenin_file=null;
            if (drawBetaCatenin)
            {
                beta_catenin_values = new double[num_lines][][]; 
                in_beta_catenin_file=new BufferedReader(new FileReader(beta_catenin_file));
                line_beta_catenin=in_beta_catenin_file.readLine();
            }
            BufferedReader in_node_file = new BufferedReader(new FileReader(node_file));
            BufferedReader in_element_file = new BufferedReader(new FileReader(element_file));
            
            String line_node = in_node_file.readLine(); // from console input example
            String line_element = in_element_file.readLine();   // above.
            
            // If line is not end of file continue
            int row = 0;
            while (line_node != null) 
            {
                // Create a StringTokenizer with a colon sign as a delimiter
                StringTokenizer st_node = new StringTokenizer(line_node);
                StringTokenizer st_element = new StringTokenizer(line_element);
                StringTokenizer st_fibre = null;
                StringTokenizer st_beta_catenin = null;
                if (drawFibres)
                {
                    st_fibre=new StringTokenizer(line_fibre);
                    Double fibre_time = Double.valueOf(st_fibre.nextToken());
                }
                
                if (drawBetaCatenin)
                {
                	st_beta_catenin=new StringTokenizer(line_beta_catenin);
                    Double beta_catenin_time = Double.valueOf(st_beta_catenin.nextToken());
                }
                Double time = Double.valueOf(st_node.nextToken());
                Double element_time = Double.valueOf(st_element.nextToken());
                if (Math.abs(time-element_time)>1e-6) 
                {
                    System.out.println("Oi - I want the element and node files with rows at the same times...");
                    System.exit(0);
                }
                
                times[row] = time.doubleValue();

                // count the number of entries in the node file and check correct 
                int entries = st_node.countTokens();
                if (entries%3 != 0)
                {
                    System.out.println("Oi - I want the node file to look like: time,x,y,type,x,y,type...");
                    System.exit(0);
                }
                numCells[row] = entries/3; 
                // count the number of entries in the element file and check correct 
                entries = st_element.countTokens();
                if (entries%3 != 0)
                {
                    System.out.println("Oi - I want the element file to look like: time,n1,n2,n3,n1,n2,n3..");
                    System.exit(0);
                }
                                
                numElements[row] = st_element.countTokens()/3;
                positions[row] = new RealPoint[memory_factor*numCells[row]];
                fibres[row] = new RealPoint[numCells[row]];
                beta_catenin_values[row] = new double[2*numCells[row]][3];
                cell_type[row]= new int[memory_factor*numCells[row]];
                element_nodes[row] = new int[memory_factor*3*numElements[row]];
                // ArrayList<Double> positionValues= new ArrayList<Double>();
                for (int i = 0; i < numCells[row]; i++) 
                {
                    double d1 = Double.valueOf(st_node.nextToken()).doubleValue();
                    double d2 = Double.valueOf(st_node.nextToken()).doubleValue();

                    if (drawFibres)
                    {
                        double f1= Double.valueOf(st_fibre.nextToken()).doubleValue();
                        double f2= Double.valueOf(st_fibre.nextToken()).doubleValue();
                        fibres[row][i]=new RealPoint(f1,f2);
                    }
                    
                    cell_type[row][i] = Integer.parseInt(st_node.nextToken());
                    if ((cell_type[row][i]<0) || (cell_type[row][i]>7))
                    {
                        System.out.println("Oi - I want a cell type between 0 and 7");
                        System.exit(0);
                    }
                    positions[row][i]=new RealPoint(d1,d2);
                }
                
                if (drawBetaCatenin)
                {
	                for (int i = 0; i < 2*numCells[row]; i++) // should not be 2*numCells[row] but due to periodic drawing don't know how many cells there will be
	                {
	                	beta_catenin_values[row][i][0]= 0.0;
	                	beta_catenin_values[row][i][1]= 0.0;
	                	beta_catenin_values[row][i][2]= 0.0;
	                }
	                //count the number of entries in the bcat file to get num non ghosts and check correct 
	                int beta_catenin_entries = st_beta_catenin.countTokens();
	                if (beta_catenin_entries%6 != 0)
	                {
	                    System.out.println("Oi - I want the node file to look like: time,index,x,y,bCat_mem,bCat_cyto,bCat_nuc,index,x,y,bCat_mem,bCat_cyto,bCat_nuc,...");
	                    System.exit(0);
	                }
	                int num_non_ghosts = beta_catenin_entries%6; 
	                for (int i = 0; i < num_non_ghosts; i++) 
	                {
                		String skip; //  Skips past unnecessary info.
                    	int index = Integer.parseInt(st_beta_catenin.nextToken()); // index
                    	skip = st_beta_catenin.nextToken(); // x
                    	skip = st_beta_catenin.nextToken(); // y
                    	
                    	double beta_catenin_membrane= Double.valueOf(st_beta_catenin.nextToken()).doubleValue();
                        double beta_catenin_cytoplasm= Double.valueOf(st_beta_catenin.nextToken()).doubleValue();
                        double beta_catenin_nuclear= Double.valueOf(st_beta_catenin.nextToken()).doubleValue();
                        beta_catenin_values[row][index][0]= beta_catenin_membrane;
                        beta_catenin_values[row][index][1]= beta_catenin_cytoplasm;
                        beta_catenin_values[row][index][2]= beta_catenin_nuclear;
                    }
                }
                for (int i = 0; i < 3*numElements[row]; i++) 
                {
                    int node = Integer.parseInt(st_element.nextToken());
                    // int node = Int.valueOf(st_element.nextToken()).intValue();
                    // positionValues.add(d1);
                    element_nodes[row][i] = node;
                }
                
                // Collections.sort(positionValues);
                // for (int i=0; i<numCells[row]; i++)
                // {
                // positions[row][i] = positionValues.get(i);
                // }

                // Read next line of the file
                if (drawFibres)
                {
                    line_fibre=in_fibre_file.readLine();
                }
                if (drawBetaCatenin)
                {
                	line_beta_catenin=in_beta_catenin_file.readLine();
                }
                line_node = in_node_file.readLine();
                line_element = in_element_file.readLine();
                row++;

            } // end while not at end of file
            
            System.out.println("Writing output files = "+writeFiles);
            System.out.println("Drawing springs = "+drawSprings);
            System.out.println("Drawing fibres = "+drawFibres);
            System.out.println("Drawing beta catenin = "+drawBetaCatenin);
            System.out.println("Drawing cells = "+drawCells);
            System.out.println("Drawing ghost nodes = "+drawGhosts);
            System.out.println("Drawing cylindrically = "+ drawCylinder);
            
            if (drawCylinder) ConvertCylindricalDataToPlane();
            
            CalculateCanvasDimensions();
            
            parsed_all_files=true;
            canvas.drawBufferedImage();
            canvas.repaint();
        } 
        catch (Exception e) 
        {
            System.out.println("Error occured. Exception message: "+e.getMessage());
        }
    }
    
    public static void CalculateCanvasDimensions()
    {
        for (int row=0 ; row<numSteps ; row++)
        {
            for (int i = 0; i < numCells[row]; i++) 
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
        for (int time_index = 0; time_index < numSteps ; time_index++)
        {
            image_cells[time_index] = new int[memory_factor*numCells[time_index]];// reserve plenty of memory
            // fill image_nodes  with an identity map (at each time step each node maps to itself)
            for (int i=0 ; i<numCells[time_index] ; i++) {
                image_cells[time_index][i] = i;
            }
            
            // draw elements first
            for (int i=0 ; i < numElements[time_index]; i++)
            {   
                // What nodes are we joining up?
                int indexA = element_nodes[time_index][3*i];
                int indexB = element_nodes[time_index][3*i+1];
                int indexC = element_nodes[time_index][3*i+2];
                // find the x-co-ords of each node
                RealPoint rA = positions[time_index][indexA];
                RealPoint rB = positions[time_index][indexB];
                RealPoint rC = positions[time_index][indexC];
                
                // identify edges that are oversized
                if ((Math.abs(rA.x - rB.x) > 0.75*crypt_width)
                    ||(Math.abs(rB.x - rC.x) > 0.75*crypt_width)
                    ||(Math.abs(rA.x - rC.x) > 0.75*crypt_width))
                {
                    MakeNewImageCell(time_index,indexA);
                    MakeNewImageCell(time_index,indexB);
                    MakeNewImageCell(time_index,indexC);
                    // break those elements into two separate elements
                    SplitElement(time_index,i);
                }
            }
        }
    }

    public static void SplitElement(int time_index,int element_index)
    {
        int indexA = element_nodes[time_index][3*element_index];
        int indexB = element_nodes[time_index][3*element_index+1];
        int indexC = element_nodes[time_index][3*element_index+2];
        // find the x-co-ords of each node
        RealPoint rA = positions[time_index][indexA];
        RealPoint rB = positions[time_index][indexB];
        RealPoint rC = positions[time_index][indexC];
        
        // Create a new element which contains the image of node A.
        element_nodes[time_index][3*numElements[time_index]] = image_cells[time_index][indexA];
        // Leave Node A in this element
        
        // if node B is far away 
        if (Math.abs(rA.x - rB.x) > 0.75*crypt_width)
        {
            // add node B to new element and add image of B to this element.
            element_nodes[time_index][3*numElements[time_index]+1] = indexB;
            element_nodes[time_index][3*element_index+1] = image_cells[time_index][indexB];
        }
        else
        {   // node B is still in this element - add its image to new element
            element_nodes[time_index][3*numElements[time_index]+1] = image_cells[time_index][indexB];
        }
        
        // if node C is far away 
        if (Math.abs(rA.x - rC.x) > 0.75*crypt_width)
        {   // add it to new element and add image of C to this element.
            element_nodes[time_index][3*numElements[time_index]+2] = indexC;
            element_nodes[time_index][3*element_index+2] = image_cells[time_index][indexC];
        }
        else
        {   // node C is still in this element - add its image to new element
            element_nodes[time_index][3*numElements[time_index]+2] = image_cells[time_index][indexC];
        }       
        numElements[time_index]++;
    }
    
    public static void MakeNewImageCell(int time_index, int node_index)
    {   // only make a new cell if one hasn't already been made
        if (image_cells[time_index][node_index]==node_index)
        {   // Make a new image of Cell A       
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
            cell_type[time_index][numCells[time_index]] = 7;
            
            // update the image record
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
        x=xs;
        y=ys;
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
        
        x=xs;
        y=ys;
    }
}

class CustomCanvas2D extends Canvas implements MouseMotionListener {
    private static final long serialVersionUID = 6997195399856046957L;

    Visualize2dCells vis;

    int width;

    int height;

    BufferedImage buffered_image=null;
    Graphics g2=null;
    
    boolean imageReady=false;   
    boolean imageDrawing=false;
    int node_radius = 5;

    Color garysSexySilver = new Color(238,238,238);
    Color garysSpringsSilver = new Color(200,200,200);
    Color ozzysDirtyGrey = new Color(80,80,80);
    Color purple = new Color(121,126,234);
    
    
    public CustomCanvas2D(Visualize2dCells v) {
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
            String filename=String.format("image%1$05d.png", vis.timeStep);
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
        int cycle=0;
        while (imageDrawing)
        {
            System.out.print("");
            if (cycle==100000)
            {
                System.out.print(".");
                cycle=0;
            }
            cycle++;
        }
        imageReady = false;
        
        int old_x = -1;
        int old_y = -1;
        int tick_length = 10;
        int num_ticks = 10;
        
        vis.time_slider.setValue(vis.timeStep);        
        vis.nutrient_colour_bar.setVisible(vis.drawNutrient);
        vis.beta_catenin_colour_bar.setVisible(vis.drawBetaCatenin);
            
        if (g2==null)
        {
            height = getHeight();
            width = getWidth();
            buffered_image=new BufferedImage(width, height,BufferedImage.TYPE_INT_RGB);
            //buffered_image=new Image(width, height);
             g2 = buffered_image.getGraphics();
        }
        
        g2.setColor(garysSexySilver);
        g2.fillRect(0,0,width,height);
        g2.setColor(Color.black);
                
        g2.drawString("Time = " + vis.times[vis.timeStep], 10,10);

        if(vis.drawCircles)
        {
	        //draw cell circle interiors
	        for (int i=0; i<vis.numCells[vis.timeStep]; i++ ) 
	        {
	            PlotPoint p=scale(vis.positions[vis.timeStep][i]);
	
	            int rx = (int) (0.5* width /(vis.max_x - vis.min_x));
	            int ry = (int) (0.5 * height /(vis.max_y - vis.min_y));
	            SetCellColour(i); 
	            if(vis.cell_type[vis.timeStep][i]!=7) // if not ghost
	            {
	            	g2.fillOval(p.x-rx, p.y-ry, 2*rx, 2*ry);
	            }
	        }        
	        
	        //draw cell circle boundaries
	        for (int i=0; i<vis.numCells[vis.timeStep]; i++ ) 
	        {
	            PlotPoint p=scale(vis.positions[vis.timeStep][i]);
	
	            int rx = (int) (0.5* width /(vis.max_x - vis.min_x));
	            int ry = (int) (0.5 * height /(vis.max_y - vis.min_y));
	            g2.setColor(Color.black);
	            if(vis.cell_type[vis.timeStep][i]!=7) // if not ghost
	            {
	            	g2.drawOval(p.x-rx, p.y-ry, 2*rx, 2*ry);
	            }
	        }        
        }
        
        double[] nutrient_conc = new double[500];    
        for  (int j=0 ; j < 500; j++)
        { 
        	nutrient_conc[j] = 0.0;
        }        
        if (vis.drawNutrient)
        {   
        	try
        	{
        		File nutrient_file = new File(vis.nutrient_file+"/nutrients_"+vis.timeStep+".dat");
        		BufferedReader in_nut_file = new BufferedReader(new FileReader(nutrient_file));

        		int row = 0;
        		String line_nut = in_nut_file.readLine(); 

        		while (line_nut != null) 
        		{
        			// Create a StringTokenizer with a colon sign as a delimiter
        			StringTokenizer st = new StringTokenizer(line_nut);
        			
        			int node_index = Integer.valueOf(st.nextToken());
        			Double x = Double.valueOf(st.nextToken());
        			Double y = Double.valueOf(st.nextToken());        			
        			nutrient_conc[node_index] = Double.valueOf(st.nextToken());
        			
        			line_nut = in_nut_file.readLine();
        			row++;                
        		}
        	}
        	catch (Exception e) 
        	{
        		System.out.println("Error occured. Exception message: "+e.getMessage());
        	}
        }
        
//        double[] membrane_beta_cat_conc = new double[500];
//        double[] cyto_beta_cat_conc = new double[500]; 
//        double[] nuc_beta_cat_conc = new double[500]; 
//        
//        for  (int j=0 ; j < 500; j++)
//        { 
//        	membrane_beta_cat_conc[j] = 0.0;
//        	cyto_beta_cat_conc[j] = 0.0;
//        	nuc_beta_cat_conc[j] = 0.0;
//        }        
//        if (vis.drawBetaCatenin)
//        {         
//        	try
//        	{
//        		File beta_catenin_file = new File(vis.beta_catenin_file+"/betacatenin_"+vis.timeStep+".dat");
//        		BufferedReader in_beta_catenin_file = new BufferedReader(new FileReader(beta_catenin_file));
//
//        		int row = 0;
//        		String line_beta_catenin = in_beta_catenin_file.readLine(); 
//
//        		while (line_beta_catenin != null) 
//        		{
//        			// Create a StringTokenizer with a colon sign as a delimiter
//        			StringTokenizer st = new StringTokenizer(line_beta_catenin);
//        			
//        			int node_index = Integer.valueOf(st.nextToken());
//        			Double x = Double.valueOf(st.nextToken());
//        			Double y = Double.valueOf(st.nextToken());        			
//        			membrane_beta_cat_conc[node_index] = Double.valueOf(st.nextToken());
//        			cyto_beta_cat_conc[node_index] = Double.valueOf(st.nextToken());
//        			nuc_beta_cat_conc[node_index] = Double.valueOf(st.nextToken());
//        			        			
//        			line_beta_catenin = in_beta_catenin_file.readLine();
//        			row++;                
//        		}
//        	}
//        	catch (Exception e) 
//        	{
//        		System.out.println("Error occured. Exception message: "+e.getMessage());
//        	}
//        }
        
        g2.setColor(Color.black);
        Shape original_clip=g2.getClip();
        // draw elements first
        for (int i=0 ; i < vis.numElements[vis.timeStep]; i++)
        {       
            // What nodes are we joining up?
        	int index[]=new int[3];
            index[0] = vis.element_nodes[vis.timeStep][3*i];
            index[1] = vis.element_nodes[vis.timeStep][3*i+1];
            index[2] = vis.element_nodes[vis.timeStep][3*i+2];
            RealPoint r1 = vis.positions[vis.timeStep][index[0]];
            RealPoint r2 = vis.positions[vis.timeStep][index[1]];
            RealPoint r3 = vis.positions[vis.timeStep][index[2]];
            
            RealPoint circumcentre=DrawCircumcentre(r1,r2,r3);
            
            PlotPoint plotcircumcentre = scale(circumcentre);
            
            // Where are they? and convert to integer pixels
            PlotPoint vertex[]=new PlotPoint[3];
            vertex[0] = scale(r1);
            vertex[1] = scale(r2);
            vertex[2] = scale(r3);

            PlotPoint midpoint[]=new PlotPoint[3];
            midpoint[2] = scale(new RealPoint(r1,r2));
            midpoint[0] = scale(new RealPoint(r2,r3));
            midpoint[1] = scale(new RealPoint(r3,r1));
                                
            g2.setColor(Color.black);
            
            if (vis.drawCells)
            {
                 int clipx[]=new int[3];
                 int clipy[]=new int[3];
                 for (int node=0;node<3;node++)
                 {
                	 clipx[node]=vertex[node].x;
                     clipy[node]=vertex[node].y;
                 }
                 Polygon clip=new Polygon(clipx,clipy,3);
                 boolean clip_me=false;
                 //Is circumcentre in the triangle?
                 //If not, then we'll clip the next bit of drawing to fit inside the triangle (ticket #432)
                 if (!clip.contains(new Point(plotcircumcentre.x, plotcircumcentre.y)))
                 {
                     clip_me=true;
                     g2.setClip(clip);
                 }
                 for (int node=0;node<3;node++)
                 {                	 
                     SetCellColour(index[node]);
                     int xs[]=new int[4];
                     int ys[]=new int[4];
                     xs[0]=plotcircumcentre.x;
                     ys[0]=plotcircumcentre.y;
                     xs[1]=midpoint[(node+1)%3].x;
                     ys[1]=midpoint[(node+1)%3].y;
                     xs[2]=vertex[node].x;
                     ys[2]=vertex[node].y;
                     xs[3]=midpoint[(node+2)%3].x;
                     ys[3]=midpoint[(node+2)%3].y;
                     g2.fillPolygon(xs,ys,4);
                 }
            
                g2.setColor(Color.black);
                ////         Plot cell boundary lines
                if( (vis.cell_type[vis.timeStep][index[0]]<7) && (vis.cell_type[vis.timeStep][index[1]]<7))
                {
                    g2.drawLine(midpoint[2].x, midpoint[2].y, plotcircumcentre.x, plotcircumcentre.y);
                }
                if( (vis.cell_type[vis.timeStep][index[1]]<7) && (vis.cell_type[vis.timeStep][index[2]]<7))
                {
                    g2.drawLine(midpoint[0].x, midpoint[0].y, plotcircumcentre.x, plotcircumcentre.y);
                }
                if( (vis.cell_type[vis.timeStep][index[2]]<7) && (vis.cell_type[vis.timeStep][index[0]]<7))
                {
                    g2.drawLine(midpoint[1].x, midpoint[1].y, plotcircumcentre.x, plotcircumcentre.y);
                }
                if (clip_me)
                {
                	g2.setClip(original_clip);
                }
            }            

            if (vis.drawNutrient)
            {                   
                int clipx[]=new int[3];
                int clipy[]=new int[3];
                for (int node=0;node<3;node++)
                {
                	clipx[node]=vertex[node].x;
                	clipy[node]=vertex[node].y;
                }
                Polygon clip=new Polygon(clipx,clipy,3);
                boolean clip_me=false;
                // Is circumcentre in the triangle?
                // If not, then we'll clip the next bit of drawing to fit inside the triangle (ticket #432)
                if (!clip.contains(new Point(plotcircumcentre.x, plotcircumcentre.y)))
                {
                	clip_me=true;
                	g2.setClip(clip);
                }
                
                for (int node=0;node<3;node++)
                {    
                	 SetCellNutrientColour(nutrient_conc[index[node]], index[node]);
                     int xs[]=new int[4];
                     int ys[]=new int[4];
                     xs[0]=plotcircumcentre.x;
                     ys[0]=plotcircumcentre.y;
                     xs[1]=midpoint[(node+1)%3].x;
                     ys[1]=midpoint[(node+1)%3].y;
                     xs[2]=vertex[node].x;
                     ys[2]=vertex[node].y;
                     xs[3]=midpoint[(node+2)%3].x;
                     ys[3]=midpoint[(node+2)%3].y;
                     g2.fillPolygon(xs,ys,4);
                }
           
                g2.setColor(Color.black);
                // Plot cell boundary lines
                if( (vis.cell_type[vis.timeStep][index[0]]<7) && (vis.cell_type[vis.timeStep][index[1]]<7))
                {
                    g2.drawLine(midpoint[2].x, midpoint[2].y, plotcircumcentre.x, plotcircumcentre.y);
                }
                if( (vis.cell_type[vis.timeStep][index[1]]<7) && (vis.cell_type[vis.timeStep][index[2]]<7))
                {
                    g2.drawLine(midpoint[0].x, midpoint[0].y, plotcircumcentre.x, plotcircumcentre.y);
                }
                if( (vis.cell_type[vis.timeStep][index[2]]<7) && (vis.cell_type[vis.timeStep][index[0]]<7))
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
                int clipx[]=new int[3];
                int clipy[]=new int[3];
                for (int node=0;node<3;node++)
                {
                	clipx[node]=vertex[node].x;
                	clipy[node]=vertex[node].y;
                }
                Polygon clip=new Polygon(clipx,clipy,3);
                boolean clip_me=false;
                // Is circumcentre in the triangle?
                // If not, then we'll clip the next bit of drawing to fit inside the triangle (ticket #432)
                if (!clip.contains(new Point(plotcircumcentre.x, plotcircumcentre.y)))
                {
                	clip_me=true;
                	g2.setClip(clip);
                }
                // plot cytoplasmic beta catenin levels
                for (int node=0;node<3;node++)
                {    
                	 SetCellCytoplasmicBetaCateninColour(vis.beta_catenin_values[vis.timeStep][index[node]][1], index[node]);
                     int xs[]=new int[4];
                     int ys[]=new int[4];
                     xs[0]=plotcircumcentre.x;
                     ys[0]=plotcircumcentre.y;
                     xs[1]=midpoint[(node+1)%3].x;
                     ys[1]=midpoint[(node+1)%3].y;
                     xs[2]=vertex[node].x;
                     ys[2]=vertex[node].y;
                     xs[3]=midpoint[(node+2)%3].x;
                     ys[3]=midpoint[(node+2)%3].y;
                     g2.fillPolygon(xs,ys,4);
                }
                // plot membrane-bound beta catenin levels
                if( (vis.cell_type[vis.timeStep][index[0]]<7) && (vis.cell_type[vis.timeStep][index[1]]<7))
                {
                	SetCellMembranBoundBetaCateninColour(vis.beta_catenin_values[vis.timeStep][index[0]][0], index[0]); 
                    g2.drawLine(midpoint[2].x, midpoint[2].y, plotcircumcentre.x, plotcircumcentre.y);
                }
                if( (vis.cell_type[vis.timeStep][index[1]]<7) && (vis.cell_type[vis.timeStep][index[2]]<7))
                {
                	SetCellMembranBoundBetaCateninColour(vis.beta_catenin_values[vis.timeStep][index[1]][0], index[1]); 
                    g2.drawLine(midpoint[0].x, midpoint[0].y, plotcircumcentre.x, plotcircumcentre.y);
                }
                if( (vis.cell_type[vis.timeStep][index[2]]<7) && (vis.cell_type[vis.timeStep][index[0]]<7))
                {
                	SetCellMembranBoundBetaCateninColour(vis.beta_catenin_values[vis.timeStep][index[2]][0], index[2]); 
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
                if( (vis.cell_type[vis.timeStep][index[0]]<7) && (vis.cell_type[vis.timeStep][index[1]]<7))
                {
                    g2.drawLine(vertex[0].x, vertex[0].y, vertex[1].x, vertex[1].y);
                }
                if( (vis.cell_type[vis.timeStep][index[1]]<7) && (vis.cell_type[vis.timeStep][index[2]]<7))
                {
                    g2.drawLine(vertex[1].x, vertex[1].y, vertex[2].x, vertex[2].y);
                }
                if( (vis.cell_type[vis.timeStep][index[2]]<7) && (vis.cell_type[vis.timeStep][index[0]]<7))
                {
                    g2.drawLine(vertex[2].x, vertex[2].y, vertex[0].x, vertex[0].y);
                }
                if(vis.drawGhosts)
                {
                    g2.setColor(garysSpringsSilver);
                    if( (vis.cell_type[vis.timeStep][index[0]]>=7) || (vis.cell_type[vis.timeStep][index[1]]>=7))
                    {
                        g2.drawLine(vertex[0].x, vertex[0].y, vertex[1].x, vertex[1].y);
                    }
                    if( (vis.cell_type[vis.timeStep][index[1]]>=7) || (vis.cell_type[vis.timeStep][index[2]]>=7))
                    {
                        g2.drawLine(vertex[1].x, vertex[1].y, vertex[2].x, vertex[2].y);
                    }
                    if( (vis.cell_type[vis.timeStep][index[2]]>=7) || (vis.cell_type[vis.timeStep][index[0]]>=7))
                    {
                        g2.drawLine(vertex[2].x, vertex[2].y, vertex[0].x, vertex[0].y);
                    }
                    g2.setColor(Color.black);
                }
            }
        }
        

        // draw nodes second so that dots are on top of lines
        double fibre_length=1.2*node_radius;
        for (int i = 0; i < vis.numCells[vis.timeStep]; i++ ) 
        {
        	PlotPoint p=scale(vis.positions[vis.timeStep][i]);

        	if(vis.drawBetaCatenin)
            { 
        		SetCellNuclearBetaCateninColour(vis.beta_catenin_values[vis.timeStep][i][2], i);
        	}
        	else
        	{
        		SetNodeColour(i);
        	}
        	if (!vis.drawNutrient)
        	{
        		g2.fillOval(p.x - node_radius, p.y - node_radius, 2 * node_radius, 2 * node_radius);
        	}
        	//old_x = p.x;
        	//old_y = p.y;
        	if (vis.drawFibres)
        	{
        		g2.setColor(Color.magenta);
        		RealPoint fibre=vis.fibres[vis.timeStep][i];
        		g2.drawLine(p.x, p.y, (int) (p.x+fibre_length*fibre.x), (int) (p.y-fibre_length*fibre.y) );

        	}
        }
        g2.setColor(Color.black);

        drawXAxis(tick_length, num_ticks);
        drawYAxis(tick_length, num_ticks);
        
        imageReady = true;
    }

    private void drawXAxis(int tick_length, int num_ticks) 
    {
        PlotPoint start=scale(vis.min_x,0);
        PlotPoint end=scale(vis.max_x,0);
        
        
        g2.drawLine(start.x, start.y, end.x, end.y);
        for (int i = 0; i <= num_ticks; i++) 
        {
            double x = vis.min_x + (i * (vis.max_x-vis.min_x)) / num_ticks;
            DecimalFormat df = new DecimalFormat("0.00");
            String x_2dp = df.format(x);
            

            //Tick lines!
            PlotPoint posn=scale(x,0);
            g2.drawLine(posn.x, posn.y, posn.x, posn.y+tick_length);
            g2.drawString(x_2dp, posn.x, posn.y + 2
                    * tick_length);
        }
    }
    
    private void drawYAxis(int tick_length, int num_ticks) 
    {
        PlotPoint start=scale(0,vis.min_y);
        PlotPoint end=scale(0,vis.max_y);
        g2.drawLine(start.x, start.y, end.x, end.y);
        
        for (int i = 0; i <= num_ticks; i++) 
        {
            double y = vis.min_y + (i * (vis.max_y-vis.min_y)) / num_ticks;
            DecimalFormat df = new DecimalFormat("0.00");
            String y_2dp = df.format(y);

            //Tick lines!
            PlotPoint posn=scale(0,y);
            g2.drawLine(posn.x-tick_length, posn.y, posn.x, posn.y);
            g2.drawString(y_2dp, posn.x - 4*tick_length, posn.y );
        

            //g2.drawString(y_2dp, scaleX(0.0) - 4* tick_length, scaleY(vis.max_y)-scaleY(y)+scaleY(0.0));
        }
    }


    
    PlotPoint scale(double x, double y)
    {
        //Map min_x to eps and max_x to width-eps (to allow a border)
        int eps=20;
        int ix = (int) ((x - vis.min_x) * (width-2*eps) /(vis.max_x - vis.min_x) +eps);
        int iy = (int) ((y - vis.min_y) * (height-2*eps) /(vis.max_y - vis.min_y) +eps);
        iy = height - iy; // This is because java is silly and has the y axis going down the screen.
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
    
    public void mouseMoved(MouseEvent e) {
        PlotPoint mouse_position = new PlotPoint(e.getX(), e.getY());
        RealPoint real_position = unscale(mouse_position);
        
        int nearest_index=-1;
        for (int i = 0; i < vis.numCells[vis.timeStep]; i++ ) 
        {
            int sq_dist=SquaredDistance(scale(vis.positions[vis.timeStep][i]),mouse_position);
            if (sq_dist < node_radius*node_radius)
            {
                nearest_index=i;
                break;
            }
        }
        if (nearest_index>=0)
        {
            RealPoint node_position=vis.positions[vis.timeStep][nearest_index];
            vis.nearest_label.setText("Node "+nearest_index+" is at "+ node_position.x + "  " + node_position.y);
        }
        else
        {
            vis.nearest_label.setText("");
        }
    }
    
    public void mouseDragged(MouseEvent e) {
        //Not used
    }
    
    int SquaredDistance(PlotPoint p0, PlotPoint p1)
    {
        int diffx=p0.x-p1.x;
        int diffy=p0.y-p1.y;
        
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
        if(vis.cell_type[vis.timeStep][index]==0)
        {
            // stem cell
            g2.setColor(Color.green);
        }
        else if (vis.cell_type[vis.timeStep][index]==1)
        {
            // transit cell
            g2.setColor(Color.orange);
        }
        else if (vis.cell_type[vis.timeStep][index]==2)
        {
            // differentiated cell
            g2.setColor(Color.red);
        }
        else if (vis.cell_type[vis.timeStep][index]==3)
        {
            // DANGER! early CANCER!
            g2.setColor(Color.gray);
        }
        else if (vis.cell_type[vis.timeStep][index]==4)
        {
            // DANGER! late CANCER!
            g2.setColor(Color.black);
        }
        else if (vis.cell_type[vis.timeStep][index]==5)
        {
            // LABELLED CELLS
            g2.setColor(Color.blue);
        }
        else if (vis.cell_type[vis.timeStep][index]==6)
        {
            // Apoptosis Cell
            g2.setColor(Color.black);
        }
        else if(vis.cell_type[vis.timeStep][index]==7)
        {
            // danger! sloughed - don't draw anything unless asked for
            if(!vis.drawGhosts)
            {
                g2.setColor(garysSexySilver);
            }
            else
            {
                g2.setColor(Color.lightGray);
            }
        }
    }
    
    
    void SetCellColour(int index)
    {
        if(vis.cell_type[vis.timeStep][index]==0)
        {
            // stem cell
            g2.setColor(Color.cyan);
        }
        else if (vis.cell_type[vis.timeStep][index]==1)
        {
            // transit cell
            g2.setColor(Color.yellow);
        }
        else if (vis.cell_type[vis.timeStep][index]==2)
        {
            // differentiated cell
            g2.setColor(Color.pink);
        }
        else if (vis.cell_type[vis.timeStep][index]==3)
        {
            // DANGER! early CANCER!
            g2.setColor(Color.lightGray);
        }
        else if (vis.cell_type[vis.timeStep][index]==4)
        {
            // DANGER! late CANCER!
            g2.setColor(Color.gray);
        }
        else if (vis.cell_type[vis.timeStep][index]==5)
        {
            // Labelled cell
            g2.setColor(purple);
        }
        else if(vis.cell_type[vis.timeStep][index]==6)
        {
            // Undergoing Apoptosis
            g2.setColor(ozzysDirtyGrey);
        }
        else if(vis.cell_type[vis.timeStep][index]==7)
        {
            // danger! sloughed - don't draw anything
            g2.setColor(garysSexySilver);
        }
    }
    
    void SetCellNutrientColour(double conc, int index)
    {
    	if(vis.cell_type[vis.timeStep][index]==6)
        {
            // Undergoing Apoptosis
            g2.setColor(ozzysDirtyGrey);
        }
        else if(vis.cell_type[vis.timeStep][index]==7)
        {
            // danger! sloughed - don't draw anything
            g2.setColor(garysSexySilver);
        }
        else
        {
            int r = (int)(121 - 80*conc); 
            if (r<0) 
            {
            	r=0;
            }
            if (r>255)
            { 
            	r=255;
            }
            int g = (int)(126 - 80*conc); 
            if (g<0) 
            {
            	g=0;
            }
            if (g>255)
            {
            	g=255;
            }
            int b = (int)(200 - 80*conc); 
            if (b<0) 
            {
            	b=0; 
            }
            if(b>255)
            {
            	b=255;
            }
            Color colour = new Color(r,g,b);
            g2.setColor(colour);
        }   	
    }
    
    void SetCellCytoplasmicBetaCateninColour(double conc, int index)
    {
    	if(vis.cell_type[vis.timeStep][index]==7)
        {
            // danger! sloughed - don't draw anything
            g2.setColor(garysSexySilver);
        }
    	else
    	{
    		int r = 100;
    		int b = 100;

    		int g = (int)(100 + 8*conc); 
    		if (g<0) 
    		{
    			g=0;
    		}
    		if (g>255)
    		{
    			g=255;
    		}

    		Color colour = new Color(r,g,b);
    		g2.setColor(colour);
    	}
    }
    
    void SetCellNuclearBetaCateninColour(double conc, int index)
    {
    	if(vis.cell_type[vis.timeStep][index]==7)
        {
            // danger! sloughed - don't draw anything
            g2.setColor(garysSexySilver);
        }
    	else
    	{
    		int r = 100;
    		int b = 100;

    		int g = (int)(100 + 8*conc); 
    		if (g<0) 
    		{
    			g=0;
    		}
    		if (g>255)
    		{
    			g=255;
    		}

    		Color colour = new Color(r,g,b);
    		g2.setColor(colour);
    	}
    }
    
    void SetCellMembranBoundBetaCateninColour(double conc, int index)
    {
    	if(vis.cell_type[vis.timeStep][index]==7)
        {
            // danger! sloughed - don't draw anything
            g2.setColor(garysSexySilver);
        }
    	else
    	{
    		int r = 100;
    		int b = 100;

    		int g = (int)(100 + 8*conc); 
    		if (g<0) 
    		{
    			g=0;
    		}
    		if (g>255)
    		{
    			g=255;
    		}

    		Color colour = new Color(r,g,b);
    		g2.setColor(colour);
    	}
    }
            
    
}
