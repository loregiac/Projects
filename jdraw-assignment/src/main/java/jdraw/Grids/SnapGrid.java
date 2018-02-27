package jdraw.Grids;

import java.awt.Point;

import jdraw.framework.DrawView;
import jdraw.framework.PointConstrainer;

public class SnapGrid implements PointConstrainer {

	private int SNAP = 15;
	private DrawView view;
	
	public SnapGrid(DrawView view){
		this.view = view;
	}
	
	@Override
	public Point constrainPoint(Point p) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getStepX(boolean right) {
		// TODO Auto-generated method stub
		return 1;
	}
	@Override
	public int getStepY(boolean down) {
		// TODO Auto-generated method stub
		return 1;
	}

	@Override
	public void activate() {
		// TODO Auto-generated method stub

	}

	@Override
	public void deactivate() {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseDown() {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseUp() {
		// TODO Auto-generated method stub

	}
	
}
