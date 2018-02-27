package jdraw.decorators;

import jdraw.framework.DrawContext;
import jdraw.framework.Figure;

public class AnimationDecorator extends DecoratorFigure {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private AnimationThread t;
	//private DrawContext context;
	
	public AnimationDecorator(Figure f, DrawContext context){
		super(f);
		//this.context = context;
		t = new AnimationThread(this,context);
		t.start();
	}
	
	public Figure unwrap(){
		t.stopThread();
		return super.getInner();
	}
}
