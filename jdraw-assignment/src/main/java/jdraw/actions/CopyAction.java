package jdraw.actions;


import java.awt.event.ActionEvent;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JMenu;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;

import jdraw.framework.DrawContext;
import jdraw.framework.Figure;
import jdraw.std.ClipBoard;
import jdraw.std.CopyCommand;

public class CopyAction extends AbstractAction implements MenuListener{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private DrawContext context;
	
	public CopyAction(DrawContext context, JMenu menu){
		this.putValue(Action.NAME,"Copy");
		//this.putValue(ACCELERATOR_KEY,KeyStroke.getKeyStroke("control C"));
		this.context = context;
		menu.addMenuListener(this);
	}
	@Override
	public void menuSelected(MenuEvent ignore) {
		setEnabled(context.getView().getSelection().size()>=1);
	}

	@Override
	public void menuDeselected(MenuEvent e) {
		setEnabled(true);
	}

	@Override
	public void menuCanceled(MenuEvent e) {
		setEnabled(true);

	}

	public void actionPerformed(ActionEvent ignore){
		List<Figure> selection =  context.getView().getSelection();
		ClipBoard.clear();
		for(Figure f : selection){
			Figure e = f.clone();
			ClipBoard.add(e);		
			context.getModel().getDrawCommandHandler().addCommand(new CopyCommand(e));
		}

	}
}

