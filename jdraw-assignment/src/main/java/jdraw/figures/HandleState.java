package jdraw.figures;

import java.awt.Cursor;
import java.awt.Point;
import java.awt.event.MouseEvent;

import jdraw.framework.DrawView;
import jdraw.framework.Figure;

public abstract class HandleState {
	
	protected final static int  SIZE = 6;
	protected final Figure owner;
	
	public HandleState(Figure owner){
		this.owner = owner;
	}
	
	public Figure getOwner(){
		return this.owner;
	}
	
	public abstract Cursor getCursor();
	
	public abstract Point getAnchor();
	
	public abstract void dragInteraction(int x, int y, MouseEvent e, DrawView v);
	
}


