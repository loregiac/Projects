package jdraw.actions;


import java.awt.event.ActionEvent;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JMenu;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;

import jdraw.framework.DrawContext;
import jdraw.framework.DrawModel;
import jdraw.framework.Figure;
import jdraw.std.ClipBoard;
import jdraw.std.PasteCommand;

public class PasteAction extends AbstractAction implements MenuListener{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private DrawContext context;
	
	public PasteAction(DrawContext context, JMenu menu){
		putValue(Action.NAME,"Paste");
		this.context = context;
		menu.addMenuListener(this);
	}
	@Override
	public void menuSelected(MenuEvent ignore) {
		setEnabled(ClipBoard.get().size()>=1);
	}

	@Override
	public void menuDeselected(MenuEvent e) {
	}

	@Override
	public void menuCanceled(MenuEvent e) {
	}

	public void actionPerformed(ActionEvent ignore){
		List<Figure> copied =  ClipBoard.get();
		DrawModel m = context.getModel();
		for(Figure f : copied){
			f.move(25, 25);
			m.addFigure(f);		
			context.getView().addToSelection(f);
			ClipBoard.remove(f);
			ClipBoard.add(f.clone());
			context.getModel().getDrawCommandHandler().addCommand(new PasteCommand(f,context));
		}
	}
}

