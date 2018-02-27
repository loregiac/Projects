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
import jdraw.std.UngroupCommand;

public class UngroupAction extends AbstractAction implements MenuListener{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private DrawContext context;
	
	public UngroupAction(DrawContext context, JMenu menu){
		putValue(Action.NAME,"Ungroup");
		this.context = context;
		menu.addMenuListener(this);
	}
	@Override
	public void menuSelected(MenuEvent ignore) {
		List<Figure> selection =  context.getView().getSelection();
		setEnabled(false);
		if(selection != null){
			for(Figure g : selection){
				if(g instanceof GroupFigure){
					setEnabled(true);
					}
			}	
		}
	}

	@Override
	public void menuDeselected(MenuEvent e) {

	}

	@Override
	public void menuCanceled(MenuEvent e) {
	}

	public void actionPerformed(ActionEvent ignore){
		List<Figure> selection =  context.getView().getSelection();
		List<Figure> totFiguresList = new LinkedList<Figure>();
		if(selection != null){
			DrawModel m = context.getView().getModel();
			for(Figure g : selection){
				if(g instanceof GroupFigure){
					Iterable<Figure> parts = ((GroupFigure) g).getFigureParts();
					m.removeFigure(g);
					context.getView().removeFromSelection(g);
					for(Figure f : parts){
						totFiguresList.add(f);
						m.addFigure(f);
						context.getView().addToSelection(f);
					}
				m.getDrawCommandHandler().addCommand(new UngroupCommand(totFiguresList,(GroupFigure) g,m,context));

				}
	
			}

		}
	}
}

