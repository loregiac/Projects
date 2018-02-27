package jdraw.decorators;

import java.awt.Rectangle;

import jdraw.framework.DrawContext;
import jdraw.framework.Figure;

public class AnimationThread extends Thread {
	private Figure f;
	private volatile boolean stop = false;
	private int dx = 1, dy = 1;
	//private DrawContext context;
	
	public AnimationThread(Figure f, DrawContext context){
		this.f = f;
		this.setDaemon(true);
		//this.context = context;
	}
	
	public void stopThread(){
		stop = true;
	}
	
	@Override
	public void run(){
		while(!stop){
			Rectangle r = f.getBounds();
			if(r.x+r.width>600 && dx>0){dx = -dx;}
			if(r.x<0 && dx <= 0){dx = -dx;}
			if(r.y+r.height>400 && dy>0){dy = -dy;}
			if(r.y<0 && dy <= 0){dy = -dy;}
			f.move(dx, dy);
			try{
				Thread.sleep(50);
			}catch(InterruptedException e){
				stop = true;
			}
			//if(!context.getView().getSelection().contains(f)){
				//stop = true;
			//}
		}
	}
}
