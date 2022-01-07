DROP TABLE IF EXISTS `login`;

CREATE TABLE `login` (
  `id` smallint unsigned NOT NULL AUTO_INCREMENT,
  `username` varchar(50) NOT NULL,
  `password` varchar(50) NOT NULL,
  `data` varchar(50) NOT NULL,
  `hora` varchar(50) NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=6 DEFAULT CHARSET=utf8;

INSERT INTO `login`(`id`,`username`,`password`,`data`,`hora`) VALUES
(DEFAULT,'teste','123','qualquer','qualquer');


DROP TABLE IF EXISTS `param_eq`;
CREATE TABLE `param_eq` (
  `cod_eq` smallint unsigned NOT NULL AUTO_INCREMENT,
  `mix` varchar(10) NOT NULL,
  `miy` varchar(10) NOT NULL,
  `rox` varchar(10) NOT NULL,
  `roy` varchar(10) NOT NULL,
  `d0` varchar(10) NOT NULL,
  `d1` varchar(10) NOT NULL,
  PRIMARY KEY (`cod_eq`)
) ENGINE=InnoDB AUTO_INCREMENT=2 DEFAULT CHARSET=utf8;

INSERT INTO `param_eq`(`cod_eq`,`mix`,`miy`,`rox`,`roy`,`d0`,`d1`) VALUES
(DEFAULT,'0.01','0.02','0.01','0.02','0.005','0.2');

DROP TABLE IF EXISTS `param_integ`;
CREATE TABLE `param_integ` (
  `cod` smallint unsigned NOT NULL AUTO_INCREMENT,
  `sitios` varchar(10) NOT NULL,
  `numeq_por_sitio` varchar(10) NOT NULL,
  `alpha` varchar(10) NOT NULL,
  `vizinhos_otimos` varchar(10) NOT NULL,
  `tinicial` varchar(10) NOT NULL,
  `tfinal` varchar(10) NOT NULL,
  `passo` varchar(10) NOT NULL,
  PRIMARY KEY (`cod`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;

INSERT INTO `param_integ`(`cod`,`sitios`,`numeq_por_sitio`,`alpha`,`vizinhos_otimos`,`tinicial`,`tfinal`,`passo`) VALUES
(DEFAULT,'101','2','1','1','0','3','1');





