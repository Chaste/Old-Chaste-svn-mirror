/*
Copyright (C) Oxford University 2008

This file is part of CHASTE.

CHASTE is free software: you can redistribute it and/or modify
it under the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

CHASTE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Lesser GNU General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CHASTE.  If not, see <http://www.gnu.org/licenses/>.
*/

import java.awt.BorderLayout;
import java.awt.Button;
import java.awt.Canvas;
import java.awt.Color;
import java.awt.Container;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Graphics2D;
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

public class Visualize  implements ActionListener, AdjustmentListener, Runnable {
	
	public Frame frame = new Frame();
	public static double[] times;
	public static int [] numCells;
	public static double[][] positions;
	public static double max_x = 0.0;
	public static int timeStep = 0;
	public static int delay=50;
	private Thread updateThread;
	CustomCanvas canvas;
	Button run;
	public static int numSteps=0; 
	
	public Visualize(){
		frame.setSize(1000, 500);
		frame.setTitle("Gavaghan's goons' visualization tools (TM)");
		
	    frame.setLayout(new BorderLayout());
		canvas = new CustomCanvas(this); 
		canvas.setSize(1000,400);
		frame.add(canvas, BorderLayout.SOUTH);
		
	    addButtons(frame);
	    frame.addWindowListener(new WindowAdapter() {
	          public void windowClosing(WindowEvent e) {
	             System.exit(0);
	             }
	          }
	      );
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
        if(pressed == "Run")
        {
        	if (timeStep==numSteps-1)
        	{
        		timeStep=0;
        	}
        	if(updateThread == null)
        	{
        		run.setLabel("Pause");
        		updateThread = new Thread(this);
        		updateThread.start();
        	}
        }
        if(pressed == "Reset")
        {
        	timeStep=0;
        	canvas.repaint();
        }
        if(pressed == "Pause")
        {
        	if(updateThread != null)
        	{
        		Thread fred = updateThread;
        		updateThread = null;
        		fred.interrupt();
        		run.setLabel("Run");
        	}
        }
    }
	
    public void adjustmentValueChanged(AdjustmentEvent e) 
    {
        delay = e.getValue();
    }    
    
    public void run()
    {
    	while(updateThread != null)
    	{
    		try
    		{
    			Thread.sleep((100-delay)*1);
    		}
    		catch(InterruptedException e)
    		{
    			return;
    		}
    		if(timeStep<numSteps-1)
    		{
    			timeStep++;
    		} else {
    			
            	if(updateThread != null)
            	{
            		Thread fred = updateThread;
            		updateThread = null;
            		fred.interrupt();
            		run.setLabel("Run");
            	}
    			
    		}
    		canvas.repaint();
    		
    	}
    }
	
	public void addButtons(Frame frame) {
		Container box=Box.createHorizontalBox();
		Button quit = new Button("Quit");
        box.add(quit);
        quit.addActionListener(this); 
        run = new Button("Run");
        box.add(run);
        run.addActionListener(this); 
        Button reset = new Button("Reset");
        box.add(reset);
        reset.addActionListener(this); 
        frame.add(box, BorderLayout.NORTH);
    
    	Container box2=Box.createHorizontalBox();
    	//box2.setLayout(new BorderLayout());
    	
        Scrollbar delay_slider=new Scrollbar(Scrollbar.HORIZONTAL, delay, 1, 1, 100);
        //delay_slider.setSize(800,25);
        delay_slider.addAdjustmentListener(this);
        Label slow=new Label("Slow");
        //slow.setSize(50,25);
        Label fast=new Label("Fast");
        //fast.setSize(50,25);
        box2.add(slow, BorderLayout.WEST);
        box2.add(delay_slider);
        box2.add(fast, BorderLayout.EAST);
        frame.add(box2);
        }
	
	public static void main(String args[]) {
		Visualize vis=new Visualize();
		
		System.out.println("Copyright Gavaghan's goons (Gary Mirams, Sarah Eastburn and Joe Pitt-Francis)");
		try{
			File file = new File(args[0]);
			BufferedReader skimFile = new BufferedReader( new FileReader( file ) );
			
			int num_lines = 0;
			while ( skimFile.readLine() != null){
				num_lines++;
			}
			
			numSteps = num_lines;
			times = new double[num_lines];
			positions = new double[num_lines][];
			numCells = new int[num_lines];
			
			BufferedReader inFile = new BufferedReader( new FileReader( file ) );
			String line = inFile.readLine();  // from console input example above.
            // If line is not end of file continue
            int row = 0;
			while ( line != null ) 
            {
                // Create a StringTokenizer with a colon sign as a delimiter
                StringTokenizer st = 
                    new StringTokenizer(line);
                
                Double time = Double.valueOf(st.nextToken());
                times[row] =  time.doubleValue();
                
                numCells[row]=st.countTokens();
                positions[row]=new double[numCells[row]];
                ArrayList<Double> positionValues= new ArrayList<Double>();
                for (int i=0; i<numCells[row]; i++)
                {
                	double d1 = Double.valueOf(st.nextToken()).doubleValue(); 
                	if (d1 > max_x)
                	{
                		max_x = d1;
                	}
                	positionValues.add(d1);
                	//positions[row][i]=d1;                
                }
                Collections.sort(positionValues);
                for (int i=0; i<numCells[row]; i++)
                {
                	positions[row][i] = positionValues.get(i);
                }
            
                // Read next line of the file
                line = inFile.readLine();
                
                row++;

            } // end while not at end of file
			
		}
		catch (Exception e){
			
		}
	}
}
class CustomCanvas extends Canvas {
	private static final long serialVersionUID = 1907161019334451833L;
	Visualize vis;
	int width;
    int height;
    Graphics2D g2;

	
	public CustomCanvas (Visualize v) {
		vis = v;
		setBackground (Color.yellow);
	}

    public void paint (Graphics g) {
	    int old_x=-1;
	    int radius=5;
	    int tick_length=10;
	    int num_ticks=10;
	    
	    g2 = (Graphics2D) g;
	    height = getHeight();
	    width = getWidth();
	    g2.drawString ("Time = "+vis.times[vis.timeStep], 10, height/4);
	    for (int i = 0;i < vis.numCells[vis.timeStep];i++)
	    {
	    	int x = convert(vis.positions[vis.timeStep][i]);
	    	g2.fillOval(x-radius,height/2-radius,2*radius,2*radius);
	    	if (i != 0)
	    	{
	    		drawSpring(old_x,  x);
	    		
	    	}
	    	old_x = x;
	    }
	    g2.drawLine(convert(0.0), 3*height/4, convert(vis.max_x), 3*height/4);
	    for (int i = 0;i<=num_ticks;i++ )
	    {
	    	double x=(i*vis.max_x)/num_ticks;
	        DecimalFormat df = new DecimalFormat("0.00");
	        String x_2dp = df.format(x);

	    	g2.drawLine(convert(x), 3*height/4, convert(x), 3*height/4+tick_length);
	    	g2.drawString (x_2dp, convert(x)-10, 3*height/4+2*tick_length);
	    }
	    g2.drawString ("x", convert(vis.max_x/2.0), 3*height/4+5*tick_length);
	}
    
    int convert(double realx)
    {
    	return (int)((double)(width)/20.0*((realx*18.0)/vis.max_x + 1));
	    
    }
    
    void drawSpring(int left, int right)
    {
    	int spring_width=5;
    	int spring_coils=40; //Divisible by 4
    	for (int i=0; i<spring_coils; i+=4)
    	{
    		//Down
    		g2.drawLine(left+(i)*(right-left)/spring_coils, height/2, left+(i+1)*(right-left)/spring_coils, height/2+spring_width);
    		//Back up
    		g2.drawLine(left+(i+1)*(right-left)/spring_coils, height/2+spring_width, left+(i+2)*(right-left)/spring_coils, height/2);
    		//Up
    		g2.drawLine(left+(i+2)*(right-left)/spring_coils, height/2, left+(i+3)*(right-left)/spring_coils, height/2-spring_width);
    		//Back down
    		g2.drawLine(left+(i+3)*(right-left)/spring_coils, height/2-spring_width, left+(i+4)*(right-left)/spring_coils, height/2);
    		
    	}
    	
    }
}