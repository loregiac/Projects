package jdraw.actions;


import java.awt.event.ActionEvent;
import java.util.LinkedList;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JMenu;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;

import jdraw.figures.GroupFigure;
import jdraw.framework.DrawContext;
import jdraw.framework.DrawModel;
import jdraw.framework.Figure;
import jdraw.std.GroupCommand;

public class GroupAction extends AbstractAction implements MenuListener{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private DrawContext context;
	
	public GroupAction(DrawContext context, JMenu menu){
		putValue(Action.NAME,"Group");
		this.context = context;
		menu.addMenuListener(this);
	}
	@Override
	public void menuSelected(MenuEvent ignore) {
		setEnabled(context.getView().getSelection().size()>1);
	}

	@Override
	public void menuDeselected(MenuEvent e) {

	}

	@Override
	public void menuCanceled(MenuEvent e) {
	}

	public void actionPerformed(ActionEvent ignore){
		List<Figure> selection =  context.getView().getSelection();
		Iterable<Figure> totFigures = context.getView().getModel().getFigures();
		List<Figure> totFiguresList = new LinkedList<Figure>();
		
		//conert totFigures Iterator to a corresponding List
		for(Figure f : totFigures){
			if(selection.contains(f)){
				totFiguresList.add(f);
			}
		}


		if(totFiguresList != null && totFiguresList.size()>=2){
			GroupFigure g = new GroupFigure(totFiguresList);
			DrawModel m = context.getView().getModel();
			m.addFigure(g);
			for(Figure f : totFiguresList){
				m.removeFigure(f);
				context.getView().removeFromSelection(f);
			}
			context.getView().addToSelection(g);
			
			m.getDrawCommandHandler().addCommand(new GroupCommand(totFiguresList,g,m,context));
		}
	}
}
