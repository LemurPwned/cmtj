from pathlib import Path

import numpy as np
import pytorch_lightning as pl
import torch
import torch.multiprocessing
import torch.nn as nn
import ujson as json
import wandb
from torch.utils.data import DataLoader, Dataset, random_split

BATCH_SIZE = 128
dr = Path('data')
KEYS = {
    'Ms1': {
        "min": 0.9,
        "max": 1.1,
        "scale": 1.0,
    },
    'Ms2': {
        "min": 1.55,
        "max": 1.77,
        "scale": 1
    },
    'Ku': {
        "min": 0.6e3,
        "max": 1e3,
        "scale": 1e-3
    },
    'J': {
        'min': -1e-3,
        'max': 1e-3,
        'scale': 1e3
    }
}


class PhysicsDataset(Dataset):
    def __init__(self, data_dir, keys, transform=None):
        super().__init__()
        self.data_dir = data_dir
        self.transform = transform
        self.gt_spectra = list(self.data_dir.glob("*.npz"))
        self.gt_parameters = json.load(
            open(self.data_dir / "all_params.json", 'r'))
        self.KEYS = keys

    def __len__(self):
        return len(self.gt_spectra)

    def __getitem__(self, index):
        with np.load(self.gt_spectra[index]) as data:
            pdata = data['spectrum']
        # normalize data
        pdata = (pdata - pdata.min()) / (pdata.max() - pdata.min())
        # pdata = torch.from_numpy(pdata).float()
        parameters = self.gt_parameters[self.gt_spectra[index].name.replace(
            ".npz", "")]
        if self.transform:
            pdata = self.transform(pdata)

        # create output tensor with normalised weights
        gt_tensor = np.asarray([(parameters[k] - self.KEYS[k]['min']) /
                                (self.KEYS[k]['max'] - self.KEYS[k]['min'])
                                for k in self.KEYS])
        return {
            "spectrum": pdata,
            "gt_tensor": gt_tensor,
            "filename": self.gt_spectra[index].name
        }


class PhysicsModel(pl.LightningModule):
    def __init__(self,
                 num_output_params: int,
                 lr: float = 1e-4,
                 loss_type: str = 'mse',
                 keep_prob: float = 0.8):
        super().__init__()
        self.save_hyperparameters()
        if loss_type == 'mse':
            self.loss = nn.MSELoss()
        elif loss_type == 'huber':
            self.loss = nn.HuberLoss()
        else:
            raise ValueError(f"Unknown loss type: {loss_type}")

        self.model = nn.Sequential(
            nn.Flatten(), nn.Linear(in_features=infeat, out_features=2048),
            nn.ReLU(), nn.Dropout(p=1 - keep_prob),
            nn.Linear(in_features=2048, out_features=1024), nn.ReLU(),
            nn.Dropout(p=1 - keep_prob),
            nn.Linear(in_features=1024, out_features=512), nn.ReLU(),
            nn.Dropout(p=1 - keep_prob),
            nn.Linear(in_features=512, out_features=128), nn.ReLU(),
            nn.Dropout(p=1 - keep_prob),
            nn.Linear(in_features=128, out_features=64), nn.ReLU(),
            nn.Dropout(p=1 - keep_prob),
            nn.Linear(in_features=64, out_features=num_output_params))

        # # L1 ImgIn shape=(?, 28, 28, 1)
        # # Conv -> (?, 28, 28, 32)
        # # Pool -> (?, 14, 14, 32)
        # self.layer1 = torch.nn.Sequential(
        #     torch.nn.Conv2d(1, 32, kernel_size=3, stride=1, padding=1),
        #     torch.nn.ReLU(), torch.nn.MaxPool2d(kernel_size=2, stride=2),
        #     torch.nn.Dropout(p=1 - keep_prob))
        # # L2 ImgIn shape=(?, 14, 14, 32)
        # # Conv      ->(?, 14, 14, 64)
        # # Pool      ->(?, 7, 7, 64)
        # self.layer2 = torch.nn.Sequential(
        #     torch.nn.Conv2d(32, 64, kernel_size=3, stride=1, padding=1),
        #     torch.nn.ReLU(), torch.nn.MaxPool2d(kernel_size=2, stride=2),
        #     torch.nn.Dropout(p=1 - keep_prob))
        # # L3 ImgIn shape=(?, 7, 7, 64)
        # # Conv ->(?, 7, 7, 128)
        # # Pool ->(?, 4, 4, 128)
        # self.layer3 = torch.nn.Sequential(
        #     torch.nn.Conv2d(64, 128, kernel_size=3, stride=1, padding=1),
        #     torch.nn.ReLU(),
        #     torch.nn.MaxPool2d(kernel_size=2, stride=2, padding=1),
        #     torch.nn.Dropout(p=1 - keep_prob))

        # # L4 FC 4x4x128 inputs -> 625 outputs
        # self.fc1 = torch.nn.Linear(4 * 4 * 128, 625, bias=True)
        # torch.nn.init.xavier_uniform(self.fc1.weight)
        # self.layer4 = torch.nn.Sequential(self.fc1, torch.nn.ReLU(),
        #                                   torch.nn.Dropout(p=1 - keep_prob))
        # # L5 Final FC 625 inputs -> 10 outputs
        # self.fc2 = torch.nn.Linear(625, num_output_params, bias=True)
        # # torch.nn.init.xavier_uniform_(self.fc2.weight)  # initialize parameters

    def forward(self, x):
        return self.model(x)
        # out = self.layer1(x)
        # out = self.layer2(out)
        # out = self.layer3(out)
        # out = out.view(out.size(0), -1)  # Flatten them for FC
        # out = self.fc1(out)
        # out = self.fc2(out)
        # return out

    def training_step(self, batch, batch_idx):
        x = batch['spectrum'].float()
        y = batch['gt_tensor'].float()
        y_hat = self.forward(x)
        loss = self.loss(y_hat, y)
        return {'loss': loss}

    def validation_step(self, batch, batch_idx):
        x = batch['spectrum'].float()
        y = batch['gt_tensor'].float()
        y_hat = self.forward(x)
        loss = self.loss(y_hat, y)
        self.log("val_loss", loss, prog_bar=True)
        return {'val_loss': loss}

    def validation_epoch_end(self, outputs):
        avg_loss = torch.stack([x['val_loss'] for x in outputs]).mean()
        self.log("avg_val_loss", avg_loss, prog_bar=True)
        return {'val_loss': avg_loss}

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=self.hparams.lr)


if __name__ == "__main__":
    # setup wandb logger
    sharing_strategy = "file_system"
    torch.multiprocessing.set_sharing_strategy(sharing_strategy)
    wandb.init(project="PIMM-NetDetection")

    pl.seed_everything(42)
    transform = None
    dataset = PhysicsDataset(dr, transform=transform, keys=KEYS)
    train_set, val_set = random_split(
        dataset=dataset,
        lengths=[0.8, 0.2],
        generator=torch.Generator().manual_seed(42))

    train_dataloader = DataLoader(train_set,
                                  batch_size=BATCH_SIZE,
                                  num_workers=0,
                                  shuffle=True)
    val_dataloader = DataLoader(val_set,
                                batch_size=BATCH_SIZE,
                                shuffle=False,
                                num_workers=4)
    model = PhysicsModel(lr=1e-5, num_output_params=len(KEYS))
    trainer = pl.Trainer(
        log_every_n_steps=5,
        devices=1,
        accelerator='mps',
        max_epochs=50,
        enable_checkpointing=pl.callbacks.ModelCheckpoint(
            'checkpoints',
            save_top_k=1,
            verbose=True,
            monitor='val_loss',
            mode='min',
        ),
        logger=pl.loggers.WandbLogger(
            name='physics_model',
            project='physics_model',
            log_model=True,
        ),
        callbacks=[
            pl.callbacks.EarlyStopping(
                monitor='val_loss',
                patience=3,
                verbose=True,
                mode='min',
            ),
        ],
    )
    trainer.fit(model, train_dataloader, val_dataloaders=val_dataloader)
