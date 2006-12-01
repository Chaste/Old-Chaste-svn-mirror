import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Canvas;
import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.image.*;
import java.awt.Image;
import java.awt.Label;
import java.awt.Scrollbar;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.StringTokenizer;
import java.lang.Math;

import javax.swing.Box;
import javax.swing.JPanel;
import javax.imageio.ImageIO;

public class Visualize2dCells implements ActionListener, AdjustmentListener, Runnable {

	public Frame frame = new Frame();

	public static double[] times;

	public static int[] numCells;
	public static int[] numElements;

	public static RealPoint[][] positions;
	public static RealPoint[][] fibres;
	public static int[][] element_nodes;
	public static int[][] cell_type;

	public static double max_x = -1e12;
	public static double max_y = -1e12;
	public static double min_x =  1e12;
	public static double min_y =  1e12;
	
	public static boolean parsed_all_files=false;
	
	public static boolean drawSprings=false;
	public static boolean drawCells=true;
	public static boolean writeFiles=false;
	public static boolean drawGhosts=false;
        public static boolean drawFibres=false; 
	public static int timeStep = 0;

	public static int delay = 50;
	
	private Thread updateThread;

	static CustomCanvas2D canvas;

	Button run;

	public static int numSteps = 0;

	public Visualize2dCells() {
		frame.setSize(800, 800);
		frame.setTitle("Gavaghan's goons' visualization tools (TM)");

		frame.setLayout(new BorderLayout());
		canvas = new CustomCanvas2D(this);
		//canvas.setSize(1000, 400);
		canvas.setPreferredSize(new Dimension(frame.getWidth(),frame.getHeight()));
		addButtons(frame);
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
			}
			if (updateThread == null) {
				run.setLabel("Pause");
				updateThread = new Thread(this);
				updateThread.start();
			}
		}
		if (pressed == "Reset") {
			timeStep = 0;
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

	public void adjustmentValueChanged(AdjustmentEvent e) {
		delay = e.getValue();
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
		//Container box = Box.createHorizontalBox();
		Button quit = new Button("Quit");
		//box.add(quit);
		quit.addActionListener(this);
		run = new Button("Run");
		//box.add(run);
		run.addActionListener(this);
		Button reset = new Button("Reset");
		//box.add(reset);
		reset.addActionListener(this);
		buttonPanel.add(quit);
		buttonPanel.add(run);
		buttonPanel.add(reset);
		//frame.add(box, BorderLayout.NORTH);

		//Container box2 = Box.createHorizontalBox();
		// box2.setLayout(new BorderLayout());
		
		JPanel scrollPanel = new JPanel();
		Scrollbar delay_slider = new Scrollbar(Scrollbar.HORIZONTAL, delay, 1,
				1, 100);
		// delay_slider.setSize(800,25);
		delay_slider.setPreferredSize(new Dimension(frame.getWidth(),20));
		delay_slider.addAdjustmentListener(this);
		Label slow = new Label("Slow");
		// slow.setSize(50,25);
		Label fast = new Label("Fast");
		// fast.setSize(50,25);
//		box2.add(slow, BorderLayout.WEST);
//		box2.add(delay_slider);
//		box2.add(fast, BorderLayout.EAST);
		//frame.add(box2);
		scrollPanel.add(slow);
		scrollPanel.add(delay_slider);
		scrollPanel.add(fast);
		
		JPanel northPanel = new JPanel(new GridLayout(2,0));
		northPanel.add(buttonPanel);
		northPanel.add(scrollPanel);
		frame.add(northPanel,BorderLayout.NORTH);
	}

	public static void main(String args[]) {
	     
		System.out.println("Copyright Gavaghan's goons (Gary Mirams, Sarah Eastburn, Pras Pathmanathan, Alex Fletcher & Joe Pitt-Francis)");
			for (int i=1; i<args.length; i++)
		{
			if (args[i].equals("output"))
			{
				writeFiles = true;
								
			}	
			if (args[i].equals("springs"))
			{
				drawSprings = true;
								
			}	
			if (args[i].equals("fibres"))
			{
				drawFibres = true;
								
			}	
			if (args[i].equals("nocells"))
			{
				drawCells = false;
			}	
			if (args[i].equals("ghosts"))
			{
				drawGhosts = true;
			}

		}
		File node_file = new File(args[0]+".viznodes");
		File element_file = new File(args[0]+".vizelements");
		if (!node_file.isFile())
		{
			System.out.println("The file "+args[0]+".viznodes doesn't exist");
			return;
		}
		if (!element_file.isFile())
		{
			System.out.println("The file "+args[0]+".vizelements doesn't exist");
			return;
		}
		
	
		
		File fibre_file=null;		
		if (drawFibres){
		    fibre_file= new File(args[0]+".vizfibres");
	            if (!fibre_file.isFile())
		    {
			System.out.println("The file "+args[0]+".vizfibres doesn't exist");
			return;
		    }
		}
		System.out.println("Writing output files = "+writeFiles);
		System.out.println("Drawing springs = "+drawSprings);
		System.out.println("Drawing fibres = "+drawFibres);
		System.out.println("Drawing cells = "+drawCells);
		System.out.println("Drawing ghost nodes = "+drawGhosts);
		
		Visualize2dCells vis = new Visualize2dCells();
		
                 try {
			BufferedReader skim_node_file = new BufferedReader(new FileReader(node_file));

			int num_lines = 0;
			while (skim_node_file.readLine() != null) {
				num_lines++;
			}

			numSteps = num_lines;
			times = new double[num_lines];
			positions = new RealPoint[num_lines][];
			cell_type = new int [num_lines][];
			numCells = new int[num_lines];
			numElements = new int[num_lines];
			element_nodes = new int[num_lines][];
			fibres =  new RealPoint[num_lines][];
			String line_fibre="";
			BufferedReader in_fibre_file=null;
			if (drawFibres)
			{
			    fibres = new RealPoint[num_lines][]; 
			    in_fibre_file=new BufferedReader(new FileReader(fibre_file));
			    line_fibre=in_fibre_file.readLine();
			}
			BufferedReader in_node_file = new BufferedReader(new FileReader(node_file));
			BufferedReader in_element_file = new BufferedReader(new FileReader(element_file));
			
			String line_node = in_node_file.readLine(); // from console input example
			String line_element = in_element_file.readLine();	// above.

			// If line is not end of file continue
			int row = 0;
			while (line_node != null) {
				// Create a StringTokenizer with a colon sign as a delimiter
				StringTokenizer st_node = new StringTokenizer(line_node);
				StringTokenizer st_element = new StringTokenizer(line_element);
				StringTokenizer st_fibre=null;
				if (drawFibres)
				{
				    st_fibre=new StringTokenizer(line_fibre);
				    Double fibre_time = Double.valueOf(st_fibre.nextToken());
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
				
				positions[row] = new RealPoint[numCells[row]];
				fibres[row] = new RealPoint[numCells[row]];
				cell_type[row]= new int[numCells[row]];
				element_nodes[row] = new int[3*numElements[row]];
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
					if ((cell_type[row][i]<0) || (cell_type[row][i]>4))
					{
						System.out.println("Oi - I want a cell type between 0 and 3");
						System.exit(0);
					}

					if (d1 > max_x) 
					{
						max_x = d1;
					}
					if (d2 > max_y) 
					{
						max_y = d2;
					} 
					if (d1 < min_x) 
					{
						min_x = d1;
					}
				    if (d2 < min_y) 
					{
						min_y = d2;
					}
					positions[row][i]=new RealPoint(d1,d2);
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
				line_node = in_node_file.readLine();
				line_element = in_element_file.readLine();
				row++;

			} // end while not at end of file
         	
			parsed_all_files=true;
            
            canvas.drawBufferedImage();
            canvas.repaint();
		} catch (Exception e) {

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

class CustomCanvas2D extends Canvas {
	private static final long serialVersionUID = 6997195399856046957L;

	Visualize2dCells vis;

	int width;

	int height;

	BufferedImage buffered_image=null;
	Graphics g2=null;
	
	boolean imageReady=false;	
	boolean imageDrawing=false;

	Color garysSexySilver = new Color(216,216,231);
	
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
		} catch (Exception e)
		{
		}
	    }
	    imageDrawing = false;
	}
	
	public void drawBufferedImage() 
	{
		while (imageDrawing)
		{
			
		}
		imageReady = false;
		
		int old_x = -1;
		int old_y = -1;
		int radius = 5;
		int tick_length = 10;
		int num_ticks = 10;
		
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
		g2.setColor(Color.black);
		
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
				////		 Plot cell boundary lines
				if( (vis.cell_type[vis.timeStep][index[0]]<4) && (vis.cell_type[vis.timeStep][index[1]]<4))
				{
					g2.drawLine(midpoint[2].x, midpoint[2].y, plotcircumcentre.x, plotcircumcentre.y);
				}
				if( (vis.cell_type[vis.timeStep][index[1]]<4) && (vis.cell_type[vis.timeStep][index[2]]<4))
				{
					g2.drawLine(midpoint[0].x, midpoint[0].y, plotcircumcentre.x, plotcircumcentre.y);
				}
				if( (vis.cell_type[vis.timeStep][index[2]]<4) && (vis.cell_type[vis.timeStep][index[0]]<4))
				{
					g2.drawLine(midpoint[1].x, midpoint[1].y, plotcircumcentre.x, plotcircumcentre.y);
				}
			}
			
			if (vis.drawSprings)
			{
				// Plot lines
				if( (vis.cell_type[vis.timeStep][index[0]]<4) && (vis.cell_type[vis.timeStep][index[1]]<4))
				{
					g2.drawLine(vertex[0].x, vertex[0].y, vertex[1].x, vertex[1].y);
				}
				if( (vis.cell_type[vis.timeStep][index[1]]<4) && (vis.cell_type[vis.timeStep][index[2]]<4))
				{
					g2.drawLine(vertex[1].x, vertex[1].y, vertex[2].x, vertex[2].y);
				}
				if( (vis.cell_type[vis.timeStep][index[2]]<4) && (vis.cell_type[vis.timeStep][index[0]]<4))
				{
					g2.drawLine(vertex[2].x, vertex[2].y, vertex[0].x, vertex[0].y);
				}
			} 
		}
		
		// draw nodes second so that dots are on top of lines
		double fibre_length=1.2*radius;
		for (int i = 0; i < vis.numCells[vis.timeStep]; i++ ) 
		{
			PlotPoint p=scale(vis.positions[vis.timeStep][i]);
			
			SetNodeColour(i);
			
			g2.fillOval(p.x - radius, p.y - radius, 2 * radius, 2 * radius);
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
			// DANGER! CANCER!
			g2.setColor(Color.black);
		}
		else if(vis.cell_type[vis.timeStep][index]==4)
		{
			// danger! sloughed - don't draw anything unless asked for
			if(!vis.drawGhosts)
            {
				Color garysSexySilver = new Color(216,216,231);
				g2.setColor(garysSexySilver);
            }
			else
			{
		        g2.setColor(Color.white);
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
			// DANGER! CANCER!
			g2.setColor(Color.gray);
		}
		else if(vis.cell_type[vis.timeStep][index]==4)
		{
			// danger! sloughed - don't draw anything
			g2.setColor(garysSexySilver);
		}
	}
			
	
}
