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

import javax.swing.Box;
import javax.swing.JPanel;

public class Visualize2D implements ActionListener, AdjustmentListener, Runnable {

	public Frame frame = new Frame();

	public static double[] times;

	public static int[] numCells;

	public static double[][] positions;

	public static double max_x = 0.0;
	
	public static double max_y = 0.0;

	public static int timeStep = 0;

	public static int delay = 50;

	private Thread updateThread;

	CustomCanvas2D canvas;

	Button run;

	public static int numSteps = 0;

	public Visualize2D() {
		frame.setSize(1000, 500);
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
		Visualize2D vis = new Visualize2D();

		System.out
				.println("Copyright Gavaghan's goons (Gary Mirams, Sarah Eastburn and Joe Pitt-Francis)");
		try {
			File file = new File(args[0]);
			BufferedReader skimFile = new BufferedReader(new FileReader(file));

			int num_lines = 0;
			while (skimFile.readLine() != null) {
				num_lines++;
			}

			numSteps = num_lines;
			times = new double[num_lines];
			positions = new double[num_lines][];
			numCells = new int[num_lines];

			BufferedReader inFile = new BufferedReader(new FileReader(file));
			String line = inFile.readLine(); // from console input example
												// above.
			// If line is not end of file continue
			int row = 0;
			while (line != null) {
				// Create a StringTokenizer with a colon sign as a delimiter
				StringTokenizer st = new StringTokenizer(line);

				Double time = Double.valueOf(st.nextToken());
				times[row] = time.doubleValue();

				numCells[row] = st.countTokens();
				positions[row] = new double[numCells[row]];
				// ArrayList<Double> positionValues= new ArrayList<Double>();
				for (int i = 0; i < numCells[row]; i+=2) 
				{
					double d1 = Double.valueOf(st.nextToken()).doubleValue();
					double d2 = Double.valueOf(st.nextToken()).doubleValue();
					if (d1 > max_x) 
					{
						max_x = d1;
					}
					if (d2 > max_y) 
					{
						max_y = d2;
					}
					// positionValues.add(d1);
					positions[row][i] = d1;
					positions[row][i+1] = d2;
				}
				// Collections.sort(positionValues);
				// for (int i=0; i<numCells[row]; i++)
				// {
				// positions[row][i] = positionValues.get(i);
				// }

				// Read next line of the file
				line = inFile.readLine();

				row++;

			} // end while not at end of file

		} catch (Exception e) {

		}
	}
}

class CustomCanvas2D extends Canvas {
	private static final long serialVersionUID = 6997195399856046957L;

	Visualize2D vis;

	int width;

	int height;

	Graphics2D g2;

	public CustomCanvas2D(Visualize2D v) {
		vis = v;
		setBackground(Color.yellow);
	}

	public void paint(Graphics g) {
		int old_x = -1;
		int old_y = -1;
		int radius = 5;
		int tick_length = 10;
		int num_ticks = 10;

		g2 = (Graphics2D) g;
		height = getHeight();
		width = getWidth();
		g2.drawString("Time = " + vis.times[vis.timeStep], 10,10);
		for (int i = 0; i < vis.numCells[vis.timeStep]; i += 2) 
		{
			int x = scaleX(vis.positions[vis.timeStep][i]);
			int y = scaleY(vis.positions[vis.timeStep][i + 1]);
			g2.fillOval(x - radius, height - y - radius, 2 * radius, 2 * radius);
			if (i != 0) 
			{
				// drawSpring(old_x, x);
				g2.drawLine(old_x, height - old_y, x, height - y);
			}
			old_x = x;
			old_y = y;
		}
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
		g2.drawString("x", scaleX(vis.max_x / 2.0), scaleY(vis.max_y)+3*tick_length);
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
		g2.drawString("y", scaleY(vis.max_y / 2.0), width-5*tick_length);
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
}