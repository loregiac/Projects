package jdraw.figures;

import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import jdraw.framework.FigureHandle;
import jdraw.states.*;

public abstract class AbstractRectangularFigure extends AbstractFigure {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Use the java.awt.Rectangle in order to save/reuse code.
	 */
	private Rectangle rectangle;
	// added myself
	private ArrayList<FigureHandle> handles;
	private int dx, dy;
	
	protected AbstractRectangularFigure(Point origin) { 
		rectangle = new Rectangle(origin);
		setHandles();
	}
	
	void setHandles(){
		handles = new ArrayList<FigureHandle>(8);
		handles.add(new Handle(rectangle.x, rectangle.y,this, new NW(this))); 
		handles.add(new Handle(rectangle.x+rectangle.width, rectangle.y,this, new NE(this))); 
		handles.add(new Handle(rectangle.x, rectangle.y+rectangle.height,this, new SW(this))); 
		handles.add(new Handle(rectangle.x+rectangle.width, rectangle.y+rectangle.height,this, new SE(this))); 
		handles.add(new Handle(rectangle.x +rectangle.width/2, rectangle.y,this, new N(this))); 
		handles.add(new Handle(rectangle.x +rectangle.width/2, rectangle.y+rectangle.height,this, new S(this))); 
		handles.add(new Handle(rectangle.x, rectangle.y + rectangle.height/2,this, new W(this))); 
		handles.add(new Handle(rectangle.x +rectangle.width, rectangle.y + rectangle.height/2,this, new E(this))); 
	}
	
	
	@Override
	public abstract void draw(Graphics g);

	@Override
	public void move(int dx, int dy) {
		if(dx!=0 || dy!=0){
			this.dx = dx;
			this.dy = dy;
			rectangle.translate(dx,dy);
			propagateFigureEvent();
		}
	}

	@Override
	public boolean contains(int x, int y) {
		return rectangle.contains(x, y);
	}

	@Override
	public void setBounds(Point origin, Point corner) {
		Rectangle original = new Rectangle(rectangle);
		rectangle.setFrameFromDiagonal(origin, corner);
		// TODO notification of change
		if(!original.equals(rectangle)){
			propagateFigureEvent();
		}
	}

	@Override
	public Rectangle getBounds() {
		return rectangle.getBounds();
	}
	
	
	public void swapHorizontal(){
		Handle NW = (Handle) handles.get(0);
		Handle NE = (Handle) handles.get(1);
		Handle SW = (Handle) handles.get(2);
		Handle SE = (Handle) handles.get(3);
		Handle W = (Handle) handles.get(6);
		Handle E = (Handle) handles.get(7);
		
		HandleState NWstate = NW.getState();
		HandleState NEstate = NE.getState();
		HandleState SWstate = SW.getState();
		HandleState SEstate = SE.getState();
		HandleState Wstate = W.getState();
		HandleState Estate = E.getState();

		NW.setState(NEstate);
		NE.setState(NWstate);
		SW.setState(SEstate);
		SE.setState(SWstate);
		W.setState(Estate);
		E.setState(Wstate);
	}
	
	
	public void swapVertical(){
		Handle NW = (Handle) handles.get(0);
		Handle NE = (Handle) handles.get(1);
		Handle SW = (Handle) handles.get(2);
		Handle SE = (Handle) handles.get(3);
		Handle N = (Handle) handles.get(4);
		Handle S = (Handle) handles.get(5);
		
		HandleState NWstate = NW.getState();
		HandleState NEstate = NE.getState();
		HandleState SWstate = SW.getState();
		HandleState SEstate = SE.getState();
		HandleState Nstate = N.getState();
		HandleState Sstate = S.getState();

		NW.setState(SWstate);
		NE.setState(SEstate);
		SW.setState(NWstate);
		SE.setState(NEstate);
		N.setState(Sstate);
		S.setState(Nstate);
	}
	
	public List<FigureHandle> getHandles(){
		if(handles == null){
			handles = new ArrayList<FigureHandle>(8);
			handles.add(new Handle(rectangle.x, rectangle.y,this, new NW(this))); 
			handles.add(new Handle(rectangle.x+rectangle.width, rectangle.y,this, new NE(this))); 
			handles.add(new Handle(rectangle.x, rectangle.y+rectangle.height,this, new SW(this))); 
			handles.add(new Handle(rectangle.x+rectangle.width, rectangle.y+rectangle.height,this, new SE(this))); 
			handles.add(new Handle(rectangle.x +rectangle.width/2, rectangle.y,this, new N(this))); 
			handles.add(new Handle(rectangle.x +rectangle.width/2, rectangle.y+rectangle.height,this, new S(this))); 
			handles.add(new Handle(rectangle.x, rectangle.y + rectangle.height/2,this, new W(this))); 
			handles.add(new Handle(rectangle.x +rectangle.width, rectangle.y + rectangle.height/2,this, new E(this))); 
		}
		return Collections.unmodifiableList(handles);
	}

	@Override
	public AbstractRectangularFigure clone(){
		AbstractRectangularFigure copy = (AbstractRectangularFigure) super.clone();
		copy.rectangle = (Rectangle) copy.rectangle.clone();
		copy.handles = null;
		return copy;
	}
	
	public int getDx(){
		return dx;
	}

	public int getDy(){
		return dy;
	}
}
