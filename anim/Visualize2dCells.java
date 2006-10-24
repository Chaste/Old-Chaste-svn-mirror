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

public class Visualize2dCells implements ActionListener, AdjustmentListener, Runnable {

	public Frame frame = new Frame();

	public static double[] times;

	public static int[] numCells;
	public static int[] numElements;

	public static RealPoint[][] positions;
	public static int[][] element_nodes;
	public static int[][] cell_type;

	public static double max_x = -1e12;
	public static double max_y = -1e12;
	public static double min_x =  1e12;
	public static double min_y =  1e12;
	
	public static boolean parsed_all_files=false;

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

	public void run() {
		while (updateThread != null) {
			try {
				Thread.sleep((100 - delay) * 1);
			} catch (InterruptedException e) {
				return;
			}
			if (timeStep < numSteps - 1) {
				timeStep++;
			} else {

				if (updateThread != null) {
					Thread fred = updateThread;
					updateThread = null;
					fred.interrupt();
					run.setLabel("Run");
				}

			}
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
		Visualize2dCells vis = new Visualize2dCells();
        
		System.out.println("Copyright Gavaghan's goons (Gary Mirams, Sarah Eastburn, Pras Pathmanathan & Joe Pitt-Francis)");
		try {
			File node_file = new File(args[0]);
			File element_file = new File(args[1]);
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
				cell_type[row]= new int[numCells[row]];
				element_nodes[row] = new int[3*numElements[row]];
				// ArrayList<Double> positionValues= new ArrayList<Double>();
				for (int i = 0; i < numCells[row]; i++) 
				{
					double d1 = Double.valueOf(st_node.nextToken()).doubleValue();
					double d2 = Double.valueOf(st_node.nextToken()).doubleValue();
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
				line_node = in_node_file.readLine();
				line_element = in_element_file.readLine();
				row++;

			} // end while not at end of file
            parsed_all_files=true;
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

	Graphics2D g2;

	Color garysSexySilver = new Color(216,216,231);
	
	public CustomCanvas2D(Visualize2dCells v) {
		vis = v;
	
		setBackground(garysSexySilver);
	}

	public void paint(Graphics g) {
		int old_x = -1;
		int old_y = -1;
		int radius = 5;
		int tick_length = 10;
		int num_ticks = 10;
		
		g2 = (Graphics2D) g;
		
        if (vis.parsed_all_files == false)
        {
          	g2.drawString("Still parsing input...", 10,10);
        	return;
        }
		
		height = getHeight();
		width = getWidth();
		g2.drawString("Time = " + vis.times[vis.timeStep], 10,10);
		g2.setPaint(Color.black);
		
		// draw elements first
		for (int i=0 ; i < vis.numElements[vis.timeStep]; i++)
		{
			// What nodes are we joining up?
			int n1 = vis.element_nodes[vis.timeStep][3*i];
			int n2 = vis.element_nodes[vis.timeStep][3*i+1];
			int n3 = vis.element_nodes[vis.timeStep][3*i+2];
			
			RealPoint r1 = vis.positions[vis.timeStep][n1];
			RealPoint r2 = vis.positions[vis.timeStep][n2];
			RealPoint r3 = vis.positions[vis.timeStep][n3];
			
			RealPoint circumcentre=DrawCircumcentre(r1,r2,r3);
			
			PlotPoint plotcircumcentre = scale(circumcentre);
			
			// Where are they? and convert to integer pixels
			PlotPoint p1 = scale(r1);
			PlotPoint p2 = scale(r2);
			PlotPoint p3 = scale(r3);
			
			PlotPoint m12 = scale(new RealPoint(r1,r2));
			PlotPoint m23 = scale(new RealPoint(r2,r3));
			PlotPoint m31 = scale(new RealPoint(r3,r1));
								
			
			g2.setPaint(garysSexySilver);
            // Plot lines
			if( (vis.cell_type[vis.timeStep][n1]<4) && (vis.cell_type[vis.timeStep][n2]<4))
			{
				g2.drawLine(p1.x, height - p1.y, p2.x, height - p2.y);
				// g2.drawLine(x1, height - y1, x2, height - y2);
			}
			if( (vis.cell_type[vis.timeStep][n2]<4) && (vis.cell_type[vis.timeStep][n3]<4))
			{
				g2.drawLine(p2.x, height - p2.y, p3.x, height - p3.y);
				// g2.drawLine(x2, height - y2, x3, height - y3);
			}
			if( (vis.cell_type[vis.timeStep][n3]<4) && (vis.cell_type[vis.timeStep][n1]<4))
			{
				g2.drawLine(p3.x, height - p3.y, p1.x, height - p1.y);
				// g2.drawLine(x3, height - y3, x1, height - y1);
			}
			
			g2.setPaint(Color.black);
			//// Plot cell boundary lines
			if( (vis.cell_type[vis.timeStep][n1]<4) && (vis.cell_type[vis.timeStep][n2]<4))
			{
				g2.drawLine(m12.x, height - m12.y, plotcircumcentre.x, height - plotcircumcentre.y);
			}
			if( (vis.cell_type[vis.timeStep][n2]<4) && (vis.cell_type[vis.timeStep][n3]<4))
			{
				g2.drawLine(m23.x, height - m23.y, plotcircumcentre.x, height - plotcircumcentre.y);
			}
			if( (vis.cell_type[vis.timeStep][n3]<4) && (vis.cell_type[vis.timeStep][n1]<4))
			{
    			g2.drawLine(m31.x, height - m31.y, plotcircumcentre.x, height - plotcircumcentre.y);
			}
    		 
			
			
		}
		
		// draw nodes second so that dots are on top of lines
		for (int i = 0; i < vis.numCells[vis.timeStep]; i++ ) 
		{
			int x = scaleX(vis.positions[vis.timeStep][i].x);
			int y = scaleY(vis.positions[vis.timeStep][i].y);
			
			if(vis.cell_type[vis.timeStep][i]==0)
			{
				// stem cell
				g2.setPaint(Color.green);
			}
			else if (vis.cell_type[vis.timeStep][i]==1)
			{
				// transit cell
				g2.setPaint(Color.orange);
			}
			else if (vis.cell_type[vis.timeStep][i]==2)
			{
				// differentiated cell
				g2.setPaint(Color.red);
			}
			else if (vis.cell_type[vis.timeStep][i]==3)
			{
				// DANGER! CANCER!
				g2.setPaint(Color.black);
			}
			else if(vis.cell_type[vis.timeStep][i]==4)
			{
				// danger! sloughed - don't draw anything
				Color garysSexySilver = new Color(216,216,231);
				g2.setPaint(garysSexySilver);
			}
				

			g2.fillOval(x - radius, height - y - radius, 2 * radius, 2 * radius);
			old_x = x;
			old_y = y;
		}
		g2.setPaint(Color.black);
		
		
//		for (int i = 0; i < vis.numCells[vis.timeStep]; i += 2) 
//		{
//			int x = scaleX(vis.positions[vis.timeStep][i]);
//			int y = scaleY(vis.positions[vis.timeStep][i + 1]);
//			g2.fillOval(x - radius, height - y - radius, 2 * radius, 2 * radius);
//			if (i != 0) 
//			{
//				// drawSpring(old_x, x);
//				g2.drawLine(old_x, height - old_y, x, height - y);
//			}
//			old_x = x;
//			old_y = y;
//		}
		drawXAxis(tick_length, num_ticks);
		drawYAxis(tick_length, num_ticks);
	}

	private void drawXAxis(int tick_length, int num_ticks) 
	{
		g2.drawLine(scaleX(0.0), scaleY(vis.max_y), scaleX(vis.max_x),
				scaleY(vis.max_y));
		for (int i = 0; i <= num_ticks; i++) 
		{
			double x = (i * vis.max_x) / num_ticks;
			DecimalFormat df = new DecimalFormat("0.00");
			String x_2dp = df.format(x);

			//Tick lines!
			g2.drawLine(scaleX(x), scaleY(vis.max_y), scaleX(x), scaleY(vis.max_y)+tick_length);
			g2.drawString(x_2dp, scaleX(x), scaleY(vis.max_y) + 2
					* tick_length);
		}
		//g2.drawString("x", scaleX(vis.max_x / 2.0), scaleY(vis.max_y)+3*tick_length);
	}
	
	private void drawYAxis(int tick_length, int num_ticks) 
	{
		g2.drawLine(scaleX(0.0), scaleY(0.0), scaleX(0.0),
				scaleY(vis.max_y));
		for (int i = 0; i <= num_ticks; i++) 
		{
			double y = (i * vis.max_y) / num_ticks;
			DecimalFormat df = new DecimalFormat("0.00");
			String y_2dp = df.format(y);

			g2.drawLine(scaleX(0.0)-tick_length, scaleY(y), scaleX(0.0), scaleY(y));
			g2.drawString(y_2dp, scaleX(0.0) - 4* tick_length, scaleY(vis.max_y)-scaleY(y)+scaleY(0.0));
		}
		//g2.drawString("y", scaleY(vis.max_y / 2.0), width-5*tick_length);
	}

	int scaleY(double realy) 
	{
		return (int) ((double) (height) / 20.0 * ((realy * 18.0) / vis.max_y + 1));

	}
	
	int scaleX(double realx) 
	{
		return (int) ((double) (width) / 20.0 * ((realx * 18.0) / vis.max_x + 1));

	}

	void drawSpring(int left, int right) {
		int spring_width = 5;
		int spring_coils = 40; // Divisible by 4
		for (int i = 0; i < spring_coils; i += 4) {
			// Down
			g2.drawLine(left + (i) * (right - left) / spring_coils, height / 2,
					left + (i + 1) * (right - left) / spring_coils, height / 2
							+ spring_width);
			// Back up
			g2.drawLine(left + (i + 1) * (right - left) / spring_coils, height
					/ 2 + spring_width, left + (i + 2) * (right - left)
					/ spring_coils, height / 2);
			// Up
			g2.drawLine(left + (i + 2) * (right - left) / spring_coils,
					height / 2, left + (i + 3) * (right - left) / spring_coils,
					height / 2 - spring_width);
			// Back down
			g2.drawLine(left + (i + 3) * (right - left) / spring_coils, height
					/ 2 - spring_width, left + (i + 4) * (right - left)
					/ spring_coils, height / 2);

		}

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
		// g2.fillRect(scaleX(x_c),height - scaleY(y_c),2,2);
		
	}
	
	PlotPoint scale(RealPoint p) 
	{
		int x = (int) ((double) (width) / 20.0 * ((p.x * 18.0) / vis.max_x + 1));
		int y = (int) ((double) (height) / 20.0 * ((p.y * 18.0) / vis.max_y + 1));
		
		return (new PlotPoint(x,y));
	}
		
	
}