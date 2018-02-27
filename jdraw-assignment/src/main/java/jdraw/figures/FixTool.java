package jdraw.figures;

import java.awt.Point;
import java.util.LinkedList;
import java.util.List;

import jdraw.framework.Figure;
import jdraw.framework.FigureEvent;
import jdraw.framework.FigureListener;

public class FixTool implements FigureListener {

	static private List<Figure> fixedFigures = new LinkedList<Figure>();
	static private List<Figure> originalFigures = new LinkedList<Figure>();
	static boolean activate = true;
	
	public void toggleFixFigure(Figure f){
		if(!fixedFigures.contains(f)){
			fixedFigures.add(f);
			originalFigures.add(f.clone());
			f.addFigureListener(this);
		}else{
			int i = fixedFigures.indexOf(f);
			fixedFigures.remove(i);
			originalFigures.remove(i);
			f.removeFigureListener(this);
		}
	}

	@Override
	public void figureChanged(FigureEvent e) {
		Figure f = e.getFigure();
		int i = fixedFigures.indexOf(f);
		if(fixedFigures.contains(f) && activate){
			activate = false;
			Point origin = new Point(originalFigures.get(i).getBounds().x,originalFigures.get(i).getBounds().y);
			Point corner = new Point(originalFigures.get(i).getBounds().x + originalFigures.get(i).getBounds().width,originalFigures.get(i).getBounds().y+originalFigures.get(i).getBounds().height);
			f.setBounds(origin, corner);
		}
		activate = true;
	}

}
