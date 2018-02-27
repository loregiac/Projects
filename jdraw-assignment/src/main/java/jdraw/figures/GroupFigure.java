package jdraw.figures;

import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

import jdraw.framework.Figure;
import jdraw.framework.FigureGroup;
import jdraw.framework.FigureHandle;

import jdraw.states.NEgroup;
import jdraw.states.NWgroup;
import jdraw.states.SEgroup;
import jdraw.states.SWgroup;


public class GroupFigure extends AbstractFigure implements FigureGroup {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private List<Figure> groupedFigures = new CopyOnWriteArrayList<>();
	ArrayList<FigureHandle> handles = new ArrayList<FigureHandle>(8);
	
	public GroupFigure(List<Figure> selection){
		groupedFigures = selection;
		this.setHandles();
	}
	
	void setHandles(){
		handles.add(new Handle(getBounds().x, getBounds().y,this, new NWgroup(this))); 
		handles.add(new Handle(getBounds().x+getBounds().width, getBounds().y,this, new NEgroup(this))); 
		handles.add(new Handle(getBounds().x, getBounds().y+getBounds().height,this, new SWgroup(this))); 
		handles.add(new Handle(getBounds().x+getBounds().width, getBounds().y+getBounds().height,this, new SEgroup(this))); 
	}
	
	//Useless
	@Override
	public void swapHorizontal() {
		Handle NW = (Handle) handles.get(0);
		Handle NE = (Handle) handles.get(1);
		Handle SW = (Handle) handles.get(2);
		Handle SE = (Handle) handles.get(3);
		
		HandleState NWstate = NW.getState();
		HandleState NEstate = NE.getState();
		HandleState SWstate = SW.getState();
		HandleState SEstate = SE.getState();

		NW.setState(NEstate);
		NE.setState(NWstate);
		SW.setState(SEstate);
		SE.setState(SWstate);
	}

	//Useless
	@Override
	public void swapVertical() {
		Handle NW = (Handle) handles.get(0);
		Handle NE = (Handle) handles.get(1);
		Handle SW = (Handle) handles.get(2);
		Handle SE = (Handle) handles.get(3);
	
		HandleState NWstate = NW.getState();
		HandleState NEstate = NE.getState();
		HandleState SWstate = SW.getState();
		HandleState SEstate = SE.getState();

		NW.setState(SWstate);
		NE.setState(SEstate);
		SW.setState(NWstate);
		SE.setState(NEstate);
	}

	@Override
	public List<Figure> getFigureParts() {
		return groupedFigures;
	}

	@Override
	public void draw(Graphics g) {
		for(Figure f : groupedFigures){
			f.draw(g);
		}
	}

	@Override
	public void move(int dx, int dy) {
		for(Figure f : groupedFigures){
			f.move(dx, dy);
		}
		propagateFigureEvent();
	}

	@Override
	public boolean contains(int x, int y) {
		return getBounds().contains(x,y);
	}

	@Override
	public void setBounds(Point origin, Point corner) {
	}

	@Override
	public Rectangle getBounds() {
		Rectangle bounds = groupedFigures.iterator().next().getBounds();
		for(Figure f : groupedFigures){
			bounds.add(f.getBounds());
		}
		return bounds;
	}

	@Override
	public List<FigureHandle> getHandles() {
		if(handles == null){
			handles = new ArrayList<FigureHandle>(4);
			handles.add(new Handle(getBounds().x, getBounds().y,this, new NWgroup(this))); 
			handles.add(new Handle(getBounds().x+getBounds().width, getBounds().y,this, new NEgroup(this))); 
			handles.add(new Handle(getBounds().x, getBounds().y+getBounds().height,this, new SWgroup(this))); 
			handles.add(new Handle(getBounds().x+getBounds().width, getBounds().y+getBounds().height,this, new SEgroup(this))); 
		}
		return handles;
	}

	@Override
	public GroupFigure clone(){
		GroupFigure copy = (GroupFigure) super.clone();
		List<Figure> partsclones = new ArrayList<Figure>(groupedFigures.size()); 
		for(Figure f : groupedFigures){
			Figure e = f.clone();
			partsclones.add(e);
		}
		copy.groupedFigures = partsclones;
		copy.handles = null;
		return copy;
	}

	@Override
	public int getDx() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getDy() {
		// TODO Auto-generated method stub
		return 0;
	}

}
