package jdraw.Grids;

import java.awt.Point;

import jdraw.framework.PointConstrainer;


public class StepGrid implements PointConstrainer {
	
	private int STEP_X;
	private int STEP_Y;
	
	public StepGrid(int STEP_X, int STEP_Y){
		this.STEP_X = STEP_X;
		this.STEP_Y = STEP_Y;
	}
	
	@Override
	public Point constrainPoint(Point p) {
		int hrem = (int) (p.getX()%STEP_X);
		int vrem = (int) (p.getY()%STEP_Y);
		int gridx = 0;
		int gridy = 0;
		if(hrem > STEP_X/2){
			gridx = (int) (p.getX()-hrem+STEP_X);
		}else if(hrem<=STEP_X/2){
			gridx = (int) (p.getX()-hrem);
		}
		if(vrem > STEP_Y/2){
			gridy = (int) (p.getY()-vrem+STEP_Y);
		}else if(vrem <= STEP_Y/2){
			gridy = (int) (p.getY()-vrem);
		}
		return new Point(gridx,gridy);
	}

	@Override
	public int getStepX(boolean right) {
		return STEP_X;		
	}

	@Override
	public int getStepY(boolean down) {
		return STEP_Y;
	}

	@Override
	public void activate() {
	}

	@Override
	public void deactivate() {
	}

	@Override
	public void mouseDown(){
		// TODO Auto-generated method stub
	}

	@Override
	public void mouseUp() {
		// TODO Auto-generated method stub

	}

}
