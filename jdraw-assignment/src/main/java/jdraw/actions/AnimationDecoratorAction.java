package jdraw.actions;


import java.awt.event.ActionEvent;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JMenu;
import javax.swing.event.MenuEvent;
import javax.swing.event.MenuListener;

import jdraw.decorators.AnimationDecorator;
import jdraw.framework.DrawContext;
import jdraw.framework.Figure;

public class AnimationDecoratorAction extends AbstractAction implements MenuListener{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private DrawContext context;
	
	public AnimationDecoratorAction(DrawContext context, JMenu menu){
		this.putValue(Action.NAME,"Animate");
		this.context = context;
		menu.addMenuListener(this);
	}
	@Override
	public void menuSelected(MenuEvent ignore) {
		setEnabled(context.getView().getSelection().size()==1);
	}

	@Override
	public void menuDeselected(MenuEvent e) {
	}

	@Override
	public void menuCanceled(MenuEvent e) {
	}

	public void actionPerformed(ActionEvent ignore){
		Figure f = context.getView().getSelection().get(0);
		context.getView().removeFromSelection(f);
		context.getModel().removeFigure(f);
		AnimationDecorator decorated = new AnimationDecorator(f,context);
		context.getModel().addFigure(decorated);
		context.getView().addToSelection(decorated);
	}
}

