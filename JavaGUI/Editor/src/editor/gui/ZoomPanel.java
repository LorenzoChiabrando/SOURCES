/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package editor.gui;

import common.Action;
import common.Util;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 *
 * @author Elvio
 */
public class ZoomPanel extends javax.swing.JPanel {
        
    public static final int zooms[] = { 
        12, 25, 50, 75, 100, 125, 150, 200, 400
    };
    public static final String zoomLabels[] = { 
        "12%", "25 %", "50 %", "75 %", "100 %", "125%", "150 %", "200 %", "400%"
    };
    
    private Action actionZoomIn, actionZoomOut, actionNormalZoom;
    private final ActionListener zoomInListener = new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent ae) {
            setZoom(zooms[jComboBox.getSelectedIndex() + 1]);
        }
    };
    private final ActionListener zoomOutListener = new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent ae) {
            setZoom(zooms[jComboBox.getSelectedIndex() - 1]);
        }
    };
    private final ActionListener normalZoomListener = new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent ae) {
            setZoom(100);
        }
    };
    
    /**
     * Creates new form ZoomPanel
     */
    public ZoomPanel() {
        initComponents();
        if (Util.isOSX()) {
            jComboBox.putClientProperty("JComboBox.isSquare", Boolean.TRUE);
        }
        
        for (String z : zoomLabels) {
            jComboBox.addItem(z);
        }
        jSlider.setMaximum(zoomLabels.length - 1);
        setZoom(100);
        jSlider.setPreferredSize(new Dimension(zooms.length*15, jSlider.getPreferredSize().height));
        jLabelZoom.setVisible(false);
    }

    public Action getActionZoomIn() { return actionZoomIn; }
    public void setActionZoomIn(Action a) {
        if (actionZoomIn != null)
            actionZoomIn.removeActionListener(zoomInListener);
        actionZoomIn = a; 
        actionZoomIn.addActionListener(zoomInListener);
    }
    
    public Action getActionZoomOut() { return actionZoomOut; }
    public void setActionZoomOut(Action a) { 
        if (actionZoomOut != null)
            actionZoomOut.removeActionListener(zoomOutListener);
        actionZoomOut = a; 
        actionZoomOut.addActionListener(zoomOutListener);
    }
    
    public Action getActionNormalZoom() { return actionNormalZoom; }
    public void setActionNormalZoom(Action a) { 
        if (actionNormalZoom != null)
            actionNormalZoom.removeActionListener(normalZoomListener);
        actionNormalZoom = a; 
        actionNormalZoom.addActionListener(normalZoomListener);
    }

    
    public final int getZoom() { 
        return zooms[ jComboBox.getSelectedIndex() ]; 
    }
    
    public final void setZoom(int nz) {
        int pos = -1;
        for (int i = 0; i < zooms.length; i++)
            if (nz == zooms[i]) {
                pos = i;
                break;
            }
        if (pos == -1)
            throw new RuntimeException("Zoom not allowed.");
        if (jComboBox.getSelectedIndex() != pos)
            jComboBox.setSelectedIndex(pos);
        
        if (jSlider.getValue() != pos)
            jSlider.setValue(pos);
        
        boolean canZoomIn = isEnabled() && (pos < zooms.length - 1);
        boolean canZoomOut = isEnabled() && (pos > 0);
        boolean canNormalZoom = isEnabled() && (nz != 100);
        jButtonZoomOut.setEnabled(isEnabled() && canZoomOut);
        jButtonZoomIn.setEnabled(isEnabled() && canZoomIn);        
        if (actionZoomIn != null)
            actionZoomIn.setEnabled(canZoomIn);
        if (actionZoomOut != null)
            actionZoomOut.setEnabled(canZoomOut);
        if (actionNormalZoom != null)
            actionNormalZoom.setEnabled(canNormalZoom);
    }
    
    public void addActionListener(ActionListener listener) {
        jComboBox.addActionListener(listener);
    }

    public void removeActionListener(ActionListener listener) {
        jComboBox.removeActionListener(listener);
    }

    @Override
    public void setEnabled(boolean enabled) {
        super.setEnabled(enabled);
        jButtonZoomIn.setEnabled(enabled);
        jButtonZoomOut.setEnabled(enabled);
        jComboBox.setEnabled(enabled);
        jSlider.setEnabled(enabled);
        jLabelZoom.setEnabled(enabled);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {
        java.awt.GridBagConstraints gridBagConstraints;

        resourceFactory = new editor.gui.ResourceFactory();
        jLabelZoom = new javax.swing.JLabel();
        jSlider = new javax.swing.JSlider();
        jComboBox = new javax.swing.JComboBox<String>();
        jButtonZoomOut = new javax.swing.JButton();
        jButtonZoomIn = new javax.swing.JButton();

        setFocusable(false);
        setLayout(new java.awt.GridBagLayout());

        jLabelZoom.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
        jLabelZoom.setText("   Zoom:");
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        add(jLabelZoom, gridBagConstraints);

        jSlider.setMinorTickSpacing(1);
        jSlider.setSnapToTicks(true);
        jSlider.setCursor(new java.awt.Cursor(java.awt.Cursor.DEFAULT_CURSOR));
        jSlider.setFocusable(false);
        jSlider.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jSliderStateChanged(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 2;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.weightx = 0.1;
        add(jSlider, gridBagConstraints);

        jComboBox.setFocusable(false);
        jComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jComboBoxActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 4;
        gridBagConstraints.gridy = 0;
        add(jComboBox, gridBagConstraints);

        jButtonZoomOut.setIcon(resourceFactory.getBulletToggleMinus16());
        jButtonZoomOut.setBorder(javax.swing.BorderFactory.createCompoundBorder(javax.swing.BorderFactory.createEmptyBorder(1, 1, 1, 1), javax.swing.BorderFactory.createEmptyBorder(3, 3, 3, 3)));
        jButtonZoomOut.setContentAreaFilled(false);
        jButtonZoomOut.setFocusable(false);
        jButtonZoomOut.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButtonZoomOutActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.insets = new java.awt.Insets(0, 6, 0, 0);
        add(jButtonZoomOut, gridBagConstraints);

        jButtonZoomIn.setIcon(resourceFactory.getBulletTogglePlus16());
        jButtonZoomIn.setBorder(javax.swing.BorderFactory.createCompoundBorder(javax.swing.BorderFactory.createEmptyBorder(1, 1, 1, 1), javax.swing.BorderFactory.createEmptyBorder(3, 3, 3, 3)));
        jButtonZoomIn.setContentAreaFilled(false);
        jButtonZoomIn.setFocusable(false);
        jButtonZoomIn.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButtonZoomInActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 3;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 6);
        add(jButtonZoomIn, gridBagConstraints);
    }// </editor-fold>//GEN-END:initComponents

    private void jSliderStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jSliderStateChanged
        setZoom(zooms[jSlider.getValue()]);
    }//GEN-LAST:event_jSliderStateChanged

    private void jComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jComboBoxActionPerformed
        setZoom(zooms[jComboBox.getSelectedIndex()]);
    }//GEN-LAST:event_jComboBoxActionPerformed

    private void jButtonZoomOutActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButtonZoomOutActionPerformed
        zoomOutListener.actionPerformed(null);
    }//GEN-LAST:event_jButtonZoomOutActionPerformed

    private void jButtonZoomInActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButtonZoomInActionPerformed
        zoomInListener.actionPerformed(null);
    }//GEN-LAST:event_jButtonZoomInActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton jButtonZoomIn;
    private javax.swing.JButton jButtonZoomOut;
    private javax.swing.JComboBox<String> jComboBox;
    private javax.swing.JLabel jLabelZoom;
    private javax.swing.JSlider jSlider;
    private editor.gui.ResourceFactory resourceFactory;
    // End of variables declaration//GEN-END:variables
}
