package jdraw.figures;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;

//import com.sun.javafx.geom.Rectangle;

import jdraw.framework.DrawView;
import jdraw.framework.Figure;
import jdraw.framework.FigureHandle;
import jdraw.std.SetBoundsCommand;

public class Handle implements FigureHandle{
	protected static final int HANDLE_SIZE = 6;
	private int x;
	private int y;
	private Figure owner;
	private HandleState state;
	private Rectangle oldBounds;
	private Rectangle newBounds;
	
	
	public Handle(int locationx, int locationy,Figure fig, HandleState state){
		this.x = locationx;
		this.y = locationy;
		this.owner = fig;
		this.state = state;
	}

	@Override
	public Figure getOwner() {
		return owner;
	}

	@Override
	public Point getLocation() {
		return new Point(this.x,this.y);
	}

	@Override
	public void draw(Graphics g) {
		g.setColor(Color.WHITE);
		g.fillRect(getState().getAnchor().x - HANDLE_SIZE/2, getState().getAnchor().y - HANDLE_SIZE/2, HANDLE_SIZE, HANDLE_SIZE);
		g.setColor(Color.BLACK);
		g.drawRect(getState().getAnchor().x - HANDLE_SIZE/2, getState().getAnchor().y - HANDLE_SIZE/2, HANDLE_SIZE, HANDLE_SIZE);
	}

	@Override
	public Cursor getCursor() {
		return state.getCursor();
	}

	@Override
	public boolean contains(int x, int y) {
		return x<getState().getAnchor().x+HANDLE_SIZE/2 & x>getState().getAnchor().x-HANDLE_SIZE/2 & y<getState().getAnchor().y+HANDLE_SIZE/2 & y>getState().getAnchor().y-HANDLE_SIZE/2;
	}

	@Override
	public void startInteraction(int x, int y, MouseEvent e, DrawView v) {
		oldBounds = this.owner.getBounds();
	}

	@Override
	public void dragInteraction(int x, int y, MouseEvent e, DrawView v) {
		getState().dragInteraction(x, y, e, v);
		v.getDrawContext().showStatusText("w: " + Integer.toString(owner.getBounds().width) + " h: " + Integer.toString(owner.getBounds().height));
	}

	@Override
	public void stopInteraction(int x, int y, MouseEvent e, DrawView v) {
		newBounds = this.owner.getBounds();
		v.getModel().getDrawCommandHandler().addCommand(new SetBoundsCommand(owner,oldBounds,newBounds));
	}

	@Override
	public HandleState getState(){
		return this.state;
	}

	@Override
	public void setState(HandleState Nstate) {
		this.state = Nstate;		
	}
	
}
